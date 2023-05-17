#' An R6 Class object, a collection of model for population of BiSBM netowrks
#'
#' @export

bisbmpop <- R6::R6Class(
  "bisbmpop",
  #
  public = list(
    n = NULL,
    A = NULL,
    M = NULL,
    mask = NULL, # 1 for NA and 0 for observed
    distribution = NULL,
    net_id = NULL,
    model_list = NULL, # A list of size Q1max * Q2max containing the best models
    discarded_model_list = NULL,# A list of size
                                # Q1max * Q2max * (nb_models - 1) containing
                                # the discarded models
    global_opts = NULL,
    fit_opts = NULL,
    separated_inits = NULL, # A nested list : Q1 init containing Q2 init
                            # with each entry containing list of size M
                            # storing the separated inits for the class
    exploration_order_list = NULL, # A list used to store the path taken
                                  # in the state space
    Z_init = NULL,
    free_density = NULL,
    free_mixture_row = NULL,
    free_mixture_col = NULL,
    ICL = NULL,
    sep_BiSBM = NULL, # A named list containing attributes for each sep BiSBM
    BICL = NULL,
    vbound = NULL,
    best_fit = NULL,
    logfactA = NULL,
    improved = NULL,
    moving_window_coordinates = NULL, # A list of size two containing the
                                      # coordinates of the bottom left and top
                                      # right points of the square
    old_moving_window_coordinates = NULL, # A list containing the previous
                                          # coordinates of the moving window




    #' @description
    #' Create a new instance of the bisbmpop object
    #'
    #' This class is generally called via the user function `estimate_colBiSBM`
    #'
    initialize = function(netlist = NULL,
                          net_id = NULL,
                          distribution = "bernoulli",
                          free_density = FALSE,
                          free_mixture_row = FALSE,
                          free_mixture_col = FALSE,
                          Z_init = NULL,
                          global_opts = list(),
                          fit_opts = list()) {
      self$A <- netlist

      # Computing the number of rows and cols
      self$n <- vector("list", 2)
      self$n[[1]] <- vapply(self$A, nrow, FUN.VALUE = .1)
      self$n[[2]] <- vapply(self$A, ncol, FUN.VALUE = .1)

      # Computing the NA mask
      self$M <- length(self$A)
        self$mask <- lapply(
          seq_along(self$A),
          function(m) {
            mask <- matrix(0, nrow(self$A[[m]]), ncol(self$A[[m]]))
            if (sum(is.na(self$A[[m]] > 0))) {
              mask[is.na(self$A[[m]])] <- 1
            }
            mask
          }
        )

      # Assigning network ids for printing
      if (is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }

      self$Z_init <- Z_init
      self$sep_BiSBM <- list(
        M = self$M,
        models = vector("list", self$M), # to store the models
        BICL = NULL, # To store the BICL vector
        Z = NULL
      )
      self$distribution <- distribution
      self$free_density <-  free_density
      self$free_mixture_row <- free_mixture_row
      self$free_mixture_col <- free_mixture_col
      self$global_opts <- list(
        Q1_min = 1L,
        Q1_max = floor(log(sum(self$n[[1]]))) + 2,
        Q2_min = 1L,
        Q2_max = floor(log(sum(self$n[[2]]))) + 2,
        spectral_init = TRUE,
        nb_models = 5L,
        depth = 1L, # By default we set a small depth
        plot_details = 1L,
        max_pass = 10L,
        verbosity = 0L,
        nb_cores = 1L,
        parallelization_vector = c(TRUE, TRUE, TRUE),
        compare_stored = TRUE)
      self$global_opts <- utils::modifyList(self$global_opts, global_opts)
      self$vbound <- matrix(
        rep(-Inf, self$global_opts$Q1_max * self$global_opts$Q2_max),
        nrow = self$global_opts$Q1_max,
        ncol = self$global_opts$Q2_max,
        byrow = TRUE
      )
      self$ICL <- matrix(
        rep(-Inf, self$global_opts$Q1_max * self$global_opts$Q2_max),
        nrow = self$global_opts$Q1_max,
        ncol = self$global_opts$Q2_max,
        byrow = TRUE
      )
      self$BICL <- matrix(
        rep(-Inf, self$global_opts$Q1_max * self$global_opts$Q2_max),
        nrow = self$global_opts$Q1_max,
        ncol = self$global_opts$Q2_max,
        byrow = TRUE
      )

      # Initializing the model_list
      self$model_list <- vector(
        "list",
        self$global_opts$Q1_max * self$global_opts$Q2_max
      )
      dim(self$model_list) <- c(
        self$global_opts$Q1_max,
        self$global_opts$Q2_max
      )

      # Initializing the exploration list
      exploration_order_list <- vector("list")

      # Initialising the discarded model_list
      self$discarded_model_list <- vector(
        "list",
        self$global_opts$Q1_max * self$global_opts$Q2_max
      )
      dim(self$discarded_model_list) <- c(
        self$global_opts$Q1_max,
        self$global_opts$Q2_max
      )

      if (self$distribution == "poisson") {
        self$logfactA <- vapply(
          seq_along(self$A),
          function(m) {
            sum(lfactorial(self$A[[m]]) * (1-self$mask[[m]]), na.rm = TRUE)
          },
          FUN.VALUE = .1)
      }
      self$fit_opts <- list(approx_pois = FALSE,
                            algo_ve = "fp",
                            minibatch = TRUE,
                            verbosity = 0)
      self$fit_opts <- utils::modifyList(self$fit_opts, fit_opts)
    },

    #' A method to perform the splitting of the clusters
    #'
    #' @param origin_model a model (fitBipartite object) to split from.
    #' @param axis a string to indicate if this is a
    #'                      column split or a row split. "row" or "col".
    #' @return best of the possible models tested
    split_clustering = function(origin_model, axis = "row") {

      # Initialize to prevent variable leaking
      row_clustering <- NULL
      col_clustering <- NULL
      split_Q <- origin_model$Q
      possible_models_size <- NULL

      # Store the clustering to keep
      if (axis == "row") {
        # If we are splitting on the rows
        row_clustering <- lapply(seq.int(self$M), function(m) {
          # We retrieve the clustering in line for
          # the origin model for the mth network
          split_clust(
            origin_model$A[[m]], # Incidence matrix
            origin_model$MAP$Z[[m]][[1]], # The row clustering
            origin_model$Q[1], # The number of current row clusters
            is_bipartite = TRUE
          )
        })

        possible_models_size <- origin_model$Q[1]

        split_Q[1] <- split_Q[1] + 1

        col_clustering <- lapply(
          seq.int(self$M),
          function(m) {
            origin_model$MAP$Z[[m]][[2]]
          }
        )
      } else if (axis == "col"){
        # If we are splitting on the columns
        col_clustering <- lapply(seq.int(self$M), function(m) {
          # We retrieve the clustering in column for
          # the origin model for the mth network
          split_clust(
            t(origin_model$A[[m]]), # Incidence matrix
            origin_model$MAP$Z[[m]][[2]], # The col clustering
            origin_model$Q[2], # The number of current col clusters
            is_bipartite = TRUE
          )
        })

        possible_models_size <- origin_model$Q[2]

        split_Q[2] <- split_Q[2] + 1

        row_clustering <- lapply(
          seq.int(self$M),
          function(m) {
            origin_model$MAP$Z[[m]][[1]]
          }
        )
      }
      # Fitting the Q splits and selecting the next model
      if (self$global_opts$parallelization_vector[2]) { 
        possible_models <- bettermc::mclapply(seq.int(possible_models_size),
          function(q) {
            # Once the row and col clustering are correctly split
            # they are merged
            if (axis == "row") {
              # If it's a row split
              q_th_Z_init <- lapply(seq.int(self$M), function(m) {
                list(row_clustering[[m]][[q]], col_clustering[[m]])
              })
            } else {
              # If it's a col split
              q_th_Z_init <- lapply(seq.int(self$M), function(m) {
                list(row_clustering[[m]], col_clustering[[m]][[q]])
              })
            }
          if (self$global_opts$compare_stored) {
            # Checking if the clustering has been encountered before
            # in the model list
            if (!is.null(self$model_list[[split_Q[1], split_Q[2]]])) {
              stored_Z1 <- unlist(
                lapply(
                  seq_along(self$model_list[[split_Q[1], split_Q[2]]]$Z),
                  function(m) self$model_list[[split_Q[1], split_Q[2]]]$Z[[m]][[1]]
                )
              )

              q_th_Z1 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[1]]
                )
              )

              ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

              stored_Z2 <- unlist(
                lapply(
                  seq_along(self$model_list[[split_Q[1], split_Q[2]]]$Z),
                  function(m) self$model_list[[split_Q[1], split_Q[2]]]$Z[[m]][[2]]
                )
              )

              q_th_Z2 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[2]]
                )
              )

              ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)

              #cat("\n1:", ARI_Z1, "| 2:", ARI_Z2)

              # Using ARI to account for label switching
              if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                if (self$global_opts$verbosity >= 4) {
                  cat(
                    "\nThe splitting memberships", q,
                    "is identical to the model stored in (", 
                    toString(split_Q), "). It will not be computed."
                  )
                }
                return(self$model_list[[split_Q[1], split_Q[2]]])
              }
            }
            # in the discarded model list
            if (!is.null(self$discarded_model_list[[split_Q[1], split_Q[2]]])) {
              for (discarded_model in self$discarded_model_list[[split_Q[1], split_Q[2]]]) {
                stored_Z1 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[1]]
                  )
                )

                q_th_Z1 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[1]]
                  )
                )

                ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

                stored_Z2 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[2]]
                  )
                )

                q_th_Z2 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[2]]
                  )
                )

                ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)
                if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                  if (self$global_opts$verbosity >= 4) {
                    cat(
                      "\nThe splitting memberships", q,
                      "is identical to a discarded model stored in (",
                      toString(split_Q), "). It will not be computed."
                    )
                  }
                  return(discarded_model)
                }
              }
            }
          }

            q_th_model <- fitBipartiteSBMPop$new(
              A = self$A,
              Q = split_Q,
              free_mixture_row = self$free_mixture_row,
              free_mixture_col = self$free_mixture_col,
              free_density = self$free_density,
              init_method = "given",
              distribution = self$distribution,
              Z = q_th_Z_init,
              net_id = self$net_id,
              fit_opts = self$fit_opts,
            )
            q_th_model$optimize()
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\n\tFitting ", q, "/", possible_models_size, "split for ",
                axis, "\n"
              )
            }

            # For free_mixture
            q_th_models <- list(q_th_model)

            if ((self$free_mixture_row | self$free_mixture_col) & self$M > 1) {
              # The levels of tolerance
              emptiness_levels <- seq(from = 0.01, to = 0.05, by = 0.01)

              if (self$global_opts$parallelization_vector[3]) {
                models_with_different_emptiness_levels <- bettermc::mclapply(
                  seq_along(emptiness_levels),
                  function(l) {
                    Cpi <- list()

                    # Compute the support for the rows
                    if (self$free_mixture_row) {
                      if (self$global_opts$verbosity >= 4) {
                        cat(
                          "\n\t\tFree mixture on the row, minimum pi threshold : ",
                          emptiness_levels[l]
                        )
                      }
                      Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                        q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                      },
                      FUN.VALUE = rep(TRUE, split_Q[1])
                      )
                      dim(Cpi[[1]]) <- c(split_Q[1], q_th_model$M)
                    } else {
                      # All networks use all clusters
                      Cpi[[1]] <- matrix(TRUE, nrow = split_Q[1], ncol = q_th_model$M)
                    }

                    # Compute the support for the columns
                    if (self$free_mixture_col) {
                      if (self$global_opts$verbosity >= 4) {
                        cat(
                          "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                          emptiness_levels[l]
                        )
                      }
                      Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                        q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                      },
                      FUN.VALUE = rep(TRUE, split_Q[2])
                      )
                      dim(Cpi[[2]]) <- c(split_Q[2], q_th_model$M)
                    } else {
                      # All networks use all clusters
                      Cpi[[2]] <- matrix(TRUE, nrow = split_Q[2], ncol = q_th_model$M)
                    }

                    # With the supports computed we can now clone the fitBipartite
                    # and fit wrt the supports.
                    q_th_model_with_supports <- q_th_model$clone()

                    # Adding the supports
                    q_th_model_with_supports$Cpi <- Cpi
                    q_th_model_with_supports$Calpha <-
                      tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                    q_th_model_with_supports$init_method <- "empty"
                    if (self$free_mixture_row | self$free_mixture_col |
                        self$global_opts$verbosity >= 4) {
                      cat(
                        "\n\t\tFitting with threshold : ",
                        emptiness_levels[l],
                        ". Threshold number ", l, "/", length(emptiness_levels)
                      )
                    }
                    q_th_model_with_supports$optimize()
                    q_th_model_with_supports
                  },
                  mc.cores = self$global_opts$nb_cores,
                  mc.allow.recursive = TRUE,
                  mc.silent = TRUE,
                  mc.retry = -1, # To prevent big crash
                  mc.progress = FALSE
                )
              } else {
                models_with_different_emptiness_levels <- lapply(
                  seq_along(emptiness_levels),
                  function(l) {
                    Cpi <- list()

                    # Compute the support for the rows
                    if (self$free_mixture_row) {
                      if (self$global_opts$verbosity >= 4) {
                        cat(
                          "\n\t\tFree mixture on the row, minimum pi threshold : ",
                          emptiness_levels[l]
                        )
                      }
                      Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                        q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                      },
                      FUN.VALUE = rep(TRUE, split_Q[1])
                      )
                      dim(Cpi[[1]]) <- c(split_Q[1], q_th_model$M)
                    } else {
                      # All networks use all clusters
                      Cpi[[1]] <- matrix(TRUE, nrow = split_Q[1], ncol = q_th_model$M)
                    }

                    # Compute the support for the columns
                    if (self$free_mixture_col) {
                      if (self$global_opts$verbosity >= 5) {
                        cat(
                          "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                          emptiness_levels[l]
                        )
                      }
                      Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                        q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                      },
                      FUN.VALUE = rep(TRUE, split_Q[2])
                      )
                      dim(Cpi[[2]]) <- c(split_Q[2], q_th_model$M)
                    } else {
                      # All networks use all clusters
                      Cpi[[2]] <- matrix(TRUE, nrow = split_Q[2], ncol = q_th_model$M)
                    }

                    # With the supports computed we can now clone the fitBipartite
                    # and fit wrt the supports.
                    q_th_model_with_supports <- q_th_model$clone()

                    # Adding the supports
                    q_th_model_with_supports$Cpi <- Cpi
                    q_th_model_with_supports$Calpha <-
                      tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                    q_th_model_with_supports$init_method <- "empty"
                    if (self$global_opts$verbosity >= 4) {
                      cat(
                        "\n\t\tFitting with threshold : ",
                        emptiness_levels[l],
                        ". Threshold number ", l, "/", length(emptiness_levels)
                      )
                    }
                    q_th_model_with_supports$optimize()
                    q_th_model_with_supports
                  }
                )
              }

              # Now the two lists are merged
              q_th_models <- append(
                q_th_models,
                models_with_different_emptiness_levels
              )
            }

            q_th_models # The list of models is returned (if no free_mixture it's a
            # one element list).
          },
          mc.cores = self$global_opts$nb_cores,
          mc.allow.recursive = TRUE,
          mc.silent = TRUE,
          mc.retry = -1, # To prevent big crash
          mc.progress = FALSE
        )
      } 
      # Fitting the Q splits and selecting the next model
      else { 
        possible_models <- lapply(
          seq.int(possible_models_size),
          function(q) {
            # Once the row and col clustering are correctly split
            # they are merged
            if (axis == "row") {
              # If it's a row split
              q_th_Z_init <- lapply(seq.int(self$M), function(m) {
                list(row_clustering[[m]][[q]], col_clustering[[m]])
              })
            } else {
              # If it's a col split
              q_th_Z_init <- lapply(seq.int(self$M), function(m) {
                list(row_clustering[[m]], col_clustering[[m]][[q]])
              })
            }

          if (self$global_opts$compare_stored) {
            # Checking if the clustering has been encountered before
            # in the model list
            if (!is.null(self$model_list[[split_Q[1], split_Q[2]]])) {
              stored_Z1 <- unlist(
                lapply(
                  seq_along(self$model_list[[split_Q[1], split_Q[2]]]$Z),
                  function(m) self$model_list[[split_Q[1], split_Q[2]]]$Z[[m]][[1]]
                )
              )

              q_th_Z1 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[1]]
                )
              )

              ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

              stored_Z2 <- unlist(
                lapply(
                  seq_along(self$model_list[[split_Q[1], split_Q[2]]]$Z),
                  function(m) self$model_list[[split_Q[1], split_Q[2]]]$Z[[m]][[2]]
                )
              )

              q_th_Z2 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[2]]
                )
              )

              ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)

              # cat("\n1:", ARI_Z1, "| 2:", ARI_Z2)

              # Using ARI to account for label switching
              if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                if (self$global_opts$verbosity >= 4) {
                  cat(
                    "\nThe splitting memberships", q,
                    "is identical to the model stored in (",
                    toString(split_Q), "). It will not be computed."
                  )
                }
                return(self$model_list[[split_Q[1], split_Q[2]]])
              }
            }
            # in the discarded model list
            if (!is.null(self$discarded_model_list[[split_Q[1], split_Q[2]]])) {
              for (discarded_model in self$discarded_model_list[[split_Q[1], split_Q[2]]]) {
                stored_Z1 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[1]]
                  )
                )

                q_th_Z1 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[1]]
                  )
                )

                ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

                stored_Z2 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[2]]
                  )
                )

                q_th_Z2 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[2]]
                  )
                )

                ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)
                if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                  if (self$global_opts$verbosity >= 4) {
                    cat(
                      "\nThe splitting memberships", q,
                      "is identical to a discarded model stored in (",
                      toString(split_Q), "). It will not be computed."
                    )
                  }
                  return(discarded_model)
                }
              }
            }
          }

            q_th_model <- fitBipartiteSBMPop$new(
              A = self$A,
              Q = split_Q,
              free_mixture_row = self$free_mixture_row,
              free_mixture_col = self$free_mixture_col,
              free_density = self$free_density,
              init_method = "given",
              distribution = self$distribution,
              Z = q_th_Z_init,
              net_id = self$net_id,
              fit_opts = self$fit_opts,
            )
            q_th_model$optimize()
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\n\tFitting ", q, "/", possible_models_size, "split for ",
                axis, "\n"
              )
            }

            # For free_mixture
            q_th_models <- list(q_th_model)

            if ((self$free_mixture_row | self$free_mixture_col) & self$M > 1) {
              # The levels of tolerance
              emptiness_levels <- seq(from = 0.01, to = 0.05, by = 0.01)
            if (self$global_opts$parallelization_vector[3]) {
              models_with_different_emptiness_levels <- bettermc::mclapply(
                seq_along(emptiness_levels),
                function(l) {
                  Cpi <- list()

                  # Compute the support for the rows
                  if (self$free_mixture_row) {
                    if (self$global_opts$verbosity >= 5) {
                      cat(
                        "\n\t\tFree mixture on the row, minimum pi threshold : ",
                        emptiness_levels[l]
                      )
                    }
                    Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                      q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                    },
                    FUN.VALUE = rep(TRUE, split_Q[1])
                    )
                    dim(Cpi[[1]]) <- c(split_Q[1], q_th_model$M)
                  } else {
                    # All networks use all clusters
                    Cpi[[1]] <- matrix(TRUE, nrow = split_Q[1], ncol = q_th_model$M)
                  }

                  # Compute the support for the columns
                  if (self$free_mixture_col) {
                    if (self$global_opts$verbosity >= 5) {
                      cat(
                        "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                        emptiness_levels[l]
                      )
                    }
                    Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                      q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                    },
                    FUN.VALUE = rep(TRUE, split_Q[2])
                    )
                    dim(Cpi[[2]]) <- c(split_Q[2], q_th_model$M)
                  } else {
                    # All networks use all clusters
                    Cpi[[2]] <- matrix(TRUE, nrow = split_Q[2], ncol = q_th_model$M)
                  }

                  # With the supports computed we can now clone the fitBipartite
                  # and fit wrt the supports.
                  q_th_model_with_supports <- q_th_model$clone()

                  # Adding the supports
                  q_th_model_with_supports$Cpi <- Cpi
                  q_th_model_with_supports$Calpha <-
                    tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                  q_th_model_with_supports$init_method <- "empty"
                  if (self$global_opts$verbosity >= 4) {
                    cat(
                      "\n\t\tFitting with threshold : ",
                      emptiness_levels[l],
                      ". Threshold number ", l, "/", length(emptiness_levels)
                    )
                  }
                  q_th_model_with_supports$optimize()
                  q_th_model_with_supports
                },
                mc.cores = self$global_opts$nb_cores,
                mc.allow.recursive = TRUE,
                mc.silent = TRUE,
                mc.retry = -1, # To prevent big crash
                mc.progress = FALSE
              )
            } else {
              models_with_different_emptiness_levels <- lapply(
                seq_along(emptiness_levels),
                function(l) {
                  Cpi <- list()

                  # Compute the support for the rows
                  if (self$free_mixture_row) {
                    if (self$global_opts$verbosity >= 5) {
                      cat(
                        "\n\t\tFree mixture on the row, minimum pi threshold : ",
                        emptiness_levels[l]
                      )
                    }
                    Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                      q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                    },
                    FUN.VALUE = rep(TRUE, split_Q[1])
                    )
                    dim(Cpi[[1]]) <- c(split_Q[1], q_th_model$M)
                  } else {
                    # All networks use all clusters
                    Cpi[[1]] <- matrix(TRUE, nrow = split_Q[1], ncol = q_th_model$M)
                  }

                  # Compute the support for the columns
                  if (self$free_mixture_col) {
                    if (self$global_opts$verbosity >= 5) {
                      cat(
                        "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                        emptiness_levels[l]
                      )
                    }
                    Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                      q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                    },
                    FUN.VALUE = rep(TRUE, split_Q[2])
                    )
                    dim(Cpi[[2]]) <- c(split_Q[2], q_th_model$M)
                  } else {
                    # All networks use all clusters
                    Cpi[[2]] <- matrix(TRUE, nrow = split_Q[2], ncol = q_th_model$M)
                  }

                  # With the supports computed we can now clone the fitBipartite
                  # and fit wrt the supports.
                  q_th_model_with_supports <- q_th_model$clone()

                  # Adding the supports
                  q_th_model_with_supports$Cpi <- Cpi
                  q_th_model_with_supports$Calpha <-
                    tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                  q_th_model_with_supports$init_method <- "empty"
                  if (self$global_opts$verbosity >= 4) {
                    cat(
                      "\n\t\tFitting with threshold : ",
                      emptiness_levels[l],
                      ". Threshold number ", l, "/", length(emptiness_levels)
                    )
                  }
                  q_th_model_with_supports$optimize()
                  q_th_model_with_supports
                }
              )
            }

              # Now the two lists are merged
              q_th_models <- append(
                q_th_models,
                models_with_different_emptiness_levels
              )
            }

            q_th_models # The list of models is returned (if no free_mixture it's a
            # one element list).
          }
        )
      }

      # If there is free_mixture it creates nestedness so we need to unlist
      possible_models <- append(list(), unlist(possible_models))

      # Now we fit all the models for the differents splits
      possible_models_BICLs <- lapply(seq_along(possible_models), function(s) {
        possible_models[[s]]$BICL
      })

      # The best in sense of BICL is
      if (self$global_opts$verbosity >= 4) {
        cat(
          "\nThe best ", axis, " split is: ",
          which.max(possible_models_BICLs)
        )
      }

      return(possible_models[[which.max(possible_models_BICLs)]])
    },

    #' A method to perform the merging of the clusters
    #'
    #' @param origin_model a model (fitBipartite object) to merge from.
    #' @param axis a string to indicate if this is a "row", "col"
    #' or "both" merge
    #' @return best of the possible models tested
    merge_clustering = function(origin_model, axis = "row") {

      # Initialize to prevent variable leaking
      row_clustering <- NULL
      col_clustering <- NULL
      merge_Q <- origin_model$Q
      possible_models_size <- 0

      # Store the clustering to keep
      switch(axis,
        row = {
          # If we are merging on the rows
          row_clustering <- lapply(
            seq.int(self$M),
            function(m) {
              merge_clust(origin_model$MAP$Z[[m]][[1]], merge_Q[1])
            }
          )

          possible_models_size <- choose(merge_Q[1], 2)

          merge_Q[1] <- merge_Q[1] - 1

          col_clustering <- lapply(
            seq.int(self$M),
            function(m) {
              origin_model$Z[[m]][[2]]
            }
          )
        },

        col = {
          # If we are splitting on the columns
          col_clustering <- lapply(
            seq.int(self$M),
            function(m) {
              merge_clust(origin_model$MAP$Z[[m]][[2]], merge_Q[2])
            }
          )

          possible_models_size <- choose(merge_Q[2], 2)

          merge_Q[2] <- merge_Q[2] - 1

          row_clustering <- lapply(
            seq.int(self$M),
            function(m) {
              origin_model$Z[[m]][[1]]
            }
          )


        },
      both = {
        # To perform both row and col merge
      })

      # Fitting the Q merges and selecting the next model
      if (self$global_opts$parallelization_vector[2]) {
      possible_models <- bettermc::mclapply(seq.int(possible_models_size),
      function(q) {
        # Once the row and col clustering are correctly merged
        # they are merged
        switch(axis,
          row = {
            # If it's a row merge
            q_th_Z_init <- lapply(seq.int(self$M), function(m) {
              list(row_clustering[[m]][[q]], col_clustering[[m]])
            })
          },
          col = {
            # If it's a col merge
            q_th_Z_init <- lapply(seq.int(self$M), function(m) {
              list(row_clustering[[m]], col_clustering[[m]][[q]])
            })
          }
        )

        if (self$global_opts$compare_stored) {
            # Checking if the clustering has been encountered before
            # in the model list
            if (!is.null(self$model_list[[merge_Q[1], merge_Q[2]]])) {
              stored_Z1 <- unlist(
                lapply(
                  seq_along(self$model_list[[merge_Q[1], merge_Q[2]]]$Z),
                  function(m) self$model_list[[merge_Q[1], merge_Q[2]]]$Z[[m]][[1]]
                )
              )

              q_th_Z1 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[1]]
                )
              )

              ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

              stored_Z2 <- unlist(
                lapply(
                  seq_along(self$model_list[[merge_Q[1], merge_Q[2]]]$Z),
                  function(m) self$model_list[[merge_Q[1], merge_Q[2]]]$Z[[m]][[2]]
                )
              )

              q_th_Z2 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[2]]
                )
              )

              ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)

              # cat("\n1:", ARI_Z1, "| 2:", ARI_Z2)

              # Using ARI to account for label switching
              if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                if (self$global_opts$verbosity >= 4) {
                  cat(
                    "\nThe merging memberships", q,
                    "is identical to the model stored in (",
                    toString(merge_Q), "). It will not be computed."
                  )
                }
                return(self$model_list[[merge_Q[1], merge_Q[2]]])
              }
            }
            # in the discarded model list
            if (!is.null(self$discarded_model_list[[merge_Q[1], merge_Q[2]]])) {
              for (discarded_model in self$discarded_model_list[[merge_Q[1], merge_Q[2]]]) {
                stored_Z1 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[1]]
                  )
                )

                q_th_Z1 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[1]]
                  )
                )

                ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

                stored_Z2 <- unlist(
                  lapply(
                    seq_along(discarded_model$Z),
                    function(m) discarded_model$Z[[m]][[2]]
                  )
                )

                q_th_Z2 <- unlist(
                  lapply(
                    seq_along(q_th_Z_init),
                    function(m) q_th_Z_init[[m]][[2]]
                  )
                )

                ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)
                if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                  if (self$global_opts$verbosity >= 4) {
                    cat(
                      "\nThe merging memberships", q,
                      "is identical to a discarded model stored in (",
                      toString(merge_Q), "). It will not be computed."
                    )
                  }
                  return(discarded_model)
                }
              }
            }
        }

        q_th_model <- fitBipartiteSBMPop$new(
          A = self$A,
          Q = merge_Q,
          free_mixture_row = self$free_mixture_row,
          free_mixture_col = self$free_mixture_col,
          free_density = self$free_density,
          init_method = "given",
          distribution = self$distribution,
          Z = q_th_Z_init,
          net_id = self$net_id,
          fit_opts = self$fit_opts,
        )
        if (self$global_opts$verbosity >= 4) {
          cat(
            "\n\tFitting ", q, "/", possible_models_size,
            "merge for", axis, "\n"
          )
        }
        q_th_model$optimize()

        # For free_mixture
        q_th_models <- list(q_th_model)

        emptiness_levels <- NULL
        if ((self$free_mixture_row | self$free_mixture_col)& self$M > 1) {
          # The levels of tolerance
          emptiness_levels <- seq(from = 0.01, to = 0.05, by = 0.01)
          if (self$global_opts$parallelization_vector[3]){
            models_with_different_emptiness_levels <- bettermc::mclapply(
              seq_along(emptiness_levels),
              function(l) {
                Cpi <- list()

                # Compute the support for the rows
                if (self$free_mixture_row) {
                  if (self$global_opts$verbosity >= 5) {
                    cat(
                      "\n\t\tFree mixture on the row, minimum pi threshold : ",
                      emptiness_levels[l]
                    )
                  }
                  Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                    q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                  },
                  FUN.VALUE = rep(TRUE, merge_Q[1])
                  )
                  dim(Cpi[[1]]) <- c(merge_Q[1], q_th_model$M)
                } else {
                  # All networks use all clusters
                  Cpi[[1]] <- matrix(TRUE, nrow = merge_Q[1], ncol = q_th_model$M)
                }

                # Compute the support for the columns
                if (self$free_mixture_col) {
                  if (self$global_opts$verbosity >= 5) {
                    cat(
                      "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                      emptiness_levels[l]
                    )
                  }
                  Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                    q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                  },
                  FUN.VALUE = rep(TRUE, merge_Q[2])
                  )
                  dim(Cpi[[2]]) <- c(merge_Q[2], q_th_model$M)
                } else {
                  # All networks use all clusters
                  Cpi[[2]] <- matrix(TRUE, nrow = merge_Q[2], ncol = q_th_model$M)
                }

                # With the supports computed we can now clone the fitBipartite
                # and fit wrt the supports.
                q_th_model_with_supports <- q_th_model$clone()

                # Adding the supports
                q_th_model_with_supports$Cpi <- Cpi
                q_th_model_with_supports$Calpha <-
                  tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                q_th_model_with_supports$init_method <- "empty"
                if (self$global_opts$verbosity >= 4) {
                  cat(
                    "\n\t\tFitting with threshold : ",
                    emptiness_levels[l],
                    ". Threshold number ", l, "/", length(emptiness_levels)
                  )
                }
                q_th_model_with_supports$optimize()
                q_th_model_with_supports
              },
              mc.cores = self$global_opts$nb_cores,
              mc.silent = TRUE,
              mc.retry = -1, # To prevent big crash
              mc.allow.recursive = TRUE)
          } else {
            models_with_different_emptiness_levels <- lapply(
              seq_along(emptiness_levels),
              function(l) {
                Cpi <- list()

                # Compute the support for the rows
                if (self$free_mixture_row) {
                  if (self$global_opts$verbosity >= 5) {
                    cat(
                      "\n\t\tFree mixture on the row, minimum pi threshold : ",
                      emptiness_levels[l]
                    )
                  }
                  Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                    q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                  },
                  FUN.VALUE = rep(TRUE, merge_Q[1])
                  )
                  dim(Cpi[[1]]) <- c(merge_Q[1], q_th_model$M)
                } else {
                  # All networks use all clusters
                  Cpi[[1]] <- matrix(TRUE, nrow = merge_Q[1], ncol = q_th_model$M)
                }

                # Compute the support for the columns
                if (self$free_mixture_col) {
                  if (self$global_opts$verbosity >= 5) {
                    cat(
                      "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                      emptiness_levels[l]
                    )
                  }
                  Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                    q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                  },
                  FUN.VALUE = rep(TRUE, merge_Q[2])
                  )
                  dim(Cpi[[2]]) <- c(merge_Q[2], q_th_model$M)
                } else {
                  # All networks use all clusters
                  Cpi[[2]] <- matrix(TRUE, nrow = merge_Q[2], ncol = q_th_model$M)
                }

                # With the supports computed we can now clone the fitBipartite
                # and fit wrt the supports.
                q_th_model_with_supports <- q_th_model$clone()

                # Adding the supports
                q_th_model_with_supports$Cpi <- Cpi
                q_th_model_with_supports$Calpha <-
                  tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
                q_th_model_with_supports$init_method <- "empty"
                if (self$global_opts$verbosity >= 4) {
                  cat(
                    "\n\t\tFitting with threshold : ",
                    emptiness_levels[l],
                    ". Threshold number ", l, "/", length(emptiness_levels)
                  )
                }
                q_th_model_with_supports$optimize()
                q_th_model_with_supports
              }
            )
          }
          # Now the two lists are merged
          q_th_models <- append(
            q_th_models,
            models_with_different_emptiness_levels
          )
        }

        q_th_models # The list of models is returned (if no free_mixture it's a
                    # one element list).
      },
        mc.cores = self$global_opts$nb_cores,
        mc.allow.recursive = TRUE,
        mc.silent = TRUE,
        mc.retry = -1, # To prevent big crash
        mc.progress = FALSE
      )
      } else {
      possible_models <- lapply(seq.int(possible_models_size),
      function(q) {
        # Once the row and col clustering are correctly merged
        # they are merged
        switch(axis,
          row = {
            # If it's a row merge
            q_th_Z_init <- lapply(seq.int(self$M), function(m) {
              list(row_clustering[[m]][[q]], col_clustering[[m]])
            })
          },
          col = {
            # If it's a col merge
            q_th_Z_init <- lapply(seq.int(self$M), function(m) {
              list(row_clustering[[m]], col_clustering[[m]][[q]])
            })
          }
        )

        if (self$global_opts$compare_stored) {
          # Checking if the clustering has been encountered before
          # in the model list
          if (!is.null(self$model_list[[merge_Q[1], merge_Q[2]]])) {
            stored_Z1 <- unlist(
              lapply(
                seq_along(self$model_list[[merge_Q[1], merge_Q[2]]]$Z),
                function(m) self$model_list[[merge_Q[1], merge_Q[2]]]$Z[[m]][[1]]
              )
            )

            q_th_Z1 <- unlist(
              lapply(
                seq_along(q_th_Z_init),
                function(m) q_th_Z_init[[m]][[1]]
              )
            )

            ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

            stored_Z2 <- unlist(
              lapply(
                seq_along(self$model_list[[merge_Q[1], merge_Q[2]]]$Z),
                function(m) self$model_list[[merge_Q[1], merge_Q[2]]]$Z[[m]][[2]]
              )
            )

            q_th_Z2 <- unlist(
              lapply(
                seq_along(q_th_Z_init),
                function(m) q_th_Z_init[[m]][[2]]
              )
            )

            ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)

            # cat("\n1:", ARI_Z1, "| 2:", ARI_Z2)

            # Using ARI to account for label switching
            if (ARI_Z1 == 1 && ARI_Z2 == 1) {
              if (self$global_opts$verbosity >= 4) {
                cat(
                  "\nThe merging memberships", q,
                  "is identical to the model stored in (", 
                  toString(merge_Q), "). It will not be computed."
                )
              }
              return(self$model_list[[merge_Q[1], merge_Q[2]]])
            }
          }
          # in the discarded model list
          if (!is.null(self$discarded_model_list[[merge_Q[1], merge_Q[2]]])) {
            for (discarded_model in self$discarded_model_list[[merge_Q[1], merge_Q[2]]]) {
              stored_Z1 <- unlist(
                lapply(
                  seq_along(discarded_model$Z),
                  function(m) discarded_model$Z[[m]][[1]]
                )
              )

              q_th_Z1 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[1]]
                )
              )

              ARI_Z1 <- aricode::ARI(stored_Z1, q_th_Z1)

              stored_Z2 <- unlist(
                lapply(
                  seq_along(discarded_model$Z),
                  function(m) discarded_model$Z[[m]][[2]]
                )
              )

              q_th_Z2 <- unlist(
                lapply(
                  seq_along(q_th_Z_init),
                  function(m) q_th_Z_init[[m]][[2]]
                )
              )

              ARI_Z2 <- aricode::ARI(stored_Z2, q_th_Z2)
              if (ARI_Z1 == 1 && ARI_Z2 == 1) {
                if (self$global_opts$verbosity >= 4){
                cat(
                  "\nThe merging memberships", q,
                  "is identical to a discarded model stored in (",
                  toString(merge_Q), "). It will not be computed."
                )
              }
                return(discarded_model)
              }
            }
          }
        }

        q_th_model <- fitBipartiteSBMPop$new(
          A = self$A,
          Q = merge_Q,
          free_mixture_row = self$free_mixture_row,
          free_mixture_col = self$free_mixture_col,
          free_density = self$free_density,
          init_method = "given",
          distribution = self$distribution,
          Z = q_th_Z_init,
          net_id = self$net_id,
          fit_opts = self$fit_opts,
        )
        if (self$global_opts$verbosity >= 4) {
          cat(
            "\n\tFitting ", q, "/", possible_models_size,
            "merge for", axis, "\n"
          )
        }
        q_th_model$optimize()

        # For free_mixture
        q_th_models <- list(q_th_model)

        emptiness_levels <- NULL
        if ((self$free_mixture_row | self$free_mixture_col)& self$M > 1) {
          # The levels of tolerance
          emptiness_levels <- seq(from = 0.01, to = 0.05, by = 0.01)

        if (self$global_opts$parallelization_vector[3]) {
          models_with_different_emptiness_levels <- bettermc::mclapply(
            seq_along(emptiness_levels),
            function(l) {
              Cpi <- list()

              # Compute the support for the rows
              if (self$free_mixture_row) {
                if (self$global_opts$verbosity >= 5) {
                  cat(
                    "\n\t\tFree mixture on the row, minimum pi threshold : ",
                    emptiness_levels[l]
                  )
                }
                Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                  q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                },
                FUN.VALUE = rep(TRUE, merge_Q[1])
                )
                dim(Cpi[[1]]) <- c(merge_Q[1], q_th_model$M)
              } else {
                # All networks use all clusters
                Cpi[[1]] <- matrix(TRUE, nrow = merge_Q[1], ncol = q_th_model$M)
              }

              # Compute the support for the columns
              if (self$free_mixture_col) {
                if (self$global_opts$verbosity >= 5) {
                  cat(
                    "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                    emptiness_levels[l]
                  )
                }
                Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                  q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                },
                FUN.VALUE = rep(TRUE, merge_Q[2])
                )
                dim(Cpi[[2]]) <- c(merge_Q[2], q_th_model$M)
              } else {
                # All networks use all clusters
                Cpi[[2]] <- matrix(TRUE, nrow = merge_Q[2], ncol = q_th_model$M)
              }

              # With the supports computed we can now clone the fitBipartite
              # and fit wrt the supports.
              q_th_model_with_supports <- q_th_model$clone()

              # Adding the supports
              q_th_model_with_supports$Cpi <- Cpi
              q_th_model_with_supports$Calpha <-
                tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
              q_th_model_with_supports$init_method <- "empty"
              if (self$global_opts$verbosity >= 4) {
                cat(
                  "\n\t\tFitting with threshold : ",
                  emptiness_levels[l],
                  ". Threshold number ", l, "/", length(emptiness_levels)
                )
              }
              q_th_model_with_supports$optimize()
              q_th_model_with_supports
            },
            mc.cores = self$global_opts$nb_cores,
            mc.silent = TRUE,
            mc.retry = -1, # To prevent big crash
            mc.allow.recursive = TRUE
          )
        } else {
          models_with_different_emptiness_levels <- lapply(
            seq_along(emptiness_levels),
            function(l) {
              Cpi <- list()

              # Compute the support for the rows
              if (self$free_mixture_row) {
                if (self$global_opts$verbosity >= 5) {
                  cat(
                    "\n\t\tFree mixture on the row, minimum pi threshold : ",
                    emptiness_levels[l]
                  )
                }
                Cpi[[1]] <- vapply(seq(q_th_model$M), function(m) {
                  q_th_model$pi[[m]][[1]] > emptiness_levels[l]
                },
                FUN.VALUE = rep(TRUE, merge_Q[1])
                )
                dim(Cpi[[1]]) <- c(merge_Q[1], q_th_model$M)
              } else {
                # All networks use all clusters
                Cpi[[1]] <- matrix(TRUE, nrow = merge_Q[1], ncol = q_th_model$M)
              }

              # Compute the support for the columns
              if (self$free_mixture_col) {
                if (self$global_opts$verbosity >= 5) {
                  cat(
                    "\n\t\tFree mixture on the columns, minimum rho threshold : ",
                    emptiness_levels[l]
                  )
                }
                Cpi[[2]] <- vapply(seq(q_th_model$M), function(m) {
                  q_th_model$pi[[m]][[2]] > emptiness_levels[l]
                },
                FUN.VALUE = rep(TRUE, merge_Q[2])
                )
                dim(Cpi[[2]]) <- c(merge_Q[2], q_th_model$M)
              } else {
                # All networks use all clusters
                Cpi[[2]] <- matrix(TRUE, nrow = merge_Q[2], ncol = q_th_model$M)
              }

              # With the supports computed we can now clone the fitBipartite
              # and fit wrt the supports.
              q_th_model_with_supports <- q_th_model$clone()

              # Adding the supports
              q_th_model_with_supports$Cpi <- Cpi
              q_th_model_with_supports$Calpha <-
                tcrossprod(Cpi[[1]], Cpi[[2]]) > 0
              q_th_model_with_supports$init_method <- "empty"
              if (self$global_opts$verbosity >= 4) {
                cat(
                  "\n\t\tFitting with threshold : ",
                  emptiness_levels[l],
                  ". Threshold number ", l, "/", length(emptiness_levels)
                )
              }
              q_th_model_with_supports$optimize()
              q_th_model_with_supports
            }
          )
        }

          # Now the two lists are merged
          q_th_models <- append(
            q_th_models,
            models_with_different_emptiness_levels
          )
        }

        q_th_models # The list of models is returned (if no free_mixture it's a
                    # one element list).
      })
      }

      # If there is free_mixture it creates nestedness so we need to unlist
      possible_models <- append(list(),unlist(possible_models))


      # Now we fit all the models for the differents splits
      possible_models_BICLs <- lapply(
        seq_along(possible_models),
        function(s) {
          possible_models[[s]]$BICL
        }
      )

      # The best in sense of BICL is)
      if (self$global_opts$verbosity >= 4) {
        cat("\nThe best ", axis, "merge is: ", which.max(possible_models_BICLs))
      }

      return(possible_models[[which.max(possible_models_BICLs)]])
    },

    #' A method to plot the state of space and its current exploration.
    #'
    #' @details the function takes no parameters and print a plot
    #' of the current state of the model_list.
    #' @return nothing
    state_space_plot = function(){
      if (self$global_opts$plot_details >= 1 & 
      !self$global_opts$parallelization_vector[1]) {

        # Creating an empty dataframe

        data_state_space <- as.data.frame(matrix(ncol=6, nrow=0))
        names(data_state_space) <- c("Q1", "Q2", "BICL", "isMaxBICL", "startingPoint", "clusteringComplete")

        for (i in seq.int(self$global_opts$Q1_max)) {
          for (j in seq.int(self$global_opts$Q2_max)) {
            if (!is.null(self$model_list[[i, j]])) {
              # If the model with Q1 and Q2 was seen
              # we add a line in the dataframe
              data_state_space[nrow(data_state_space) + 1, ] <- list(
                i,
                j,
                self$model_list[[i, j]]$BICL,
                FALSE,
                as.character(toString(c(self$model_list[[i, j]]$greedy_exploration_starting_point[1], self$model_list[[i, j]]$greedy_exploration_starting_point[2]))),
                self$model_list[[i,j]]$clustering_is_complete
              )
            }
          }
        }

        # data_state_space[nrow(data_state_space) + 1,] <- list(
        #   0,
        #   0,
        #   -Inf,
        #   FALSE,
        #   "",
        #   FALSE
        # )

        # Here the max BICL is highlighted
        data_state_space[which.max(data_state_space$BICL), ]$isMaxBICL <- TRUE

        # Here the greedy exploration path is computed
        exploration_path <- as.data.frame(matrix(nrow = 0, ncol = 3))
        names(exploration_path) <- c("Q1", "Q2", "startingPoint")

        for (path in seq_along(self$exploration_order_list)) {
          for (idx in seq_along(self$exploration_order_list[[path]])) {
            exploration_path[nrow(exploration_path) + 1, ] <- list(
              self$exploration_order_list[[path]][[idx]][1], # Q1
              self$exploration_order_list[[path]][[idx]][2], # Q2
              as.character(toString(c( # startingPoint
                self$exploration_order_list[[path]][[1]][1],
                self$exploration_order_list[[path]][[1]][2]
              )))
            )
          }
        }

        # Plotting
        state_plot <- ggplot(data_state_space) +
          geom_path(
            data = exploration_path,
            position = position_dodge2(width = 0.2),
            aes(
              x = Q1,
              y = Q2,
              color = startingPoint
            ),
            linewidth = 2,
            arrow = arrow()
          ) +
          guides(color = guide_legend(title = "Path taken")) +
          scale_color_discrete() +
          scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))), limits = c(1, self$global_opts$Q2_max)) +
          scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))), limits = c(1, self$global_opts$Q1_max))

          if (!is.null(self$old_moving_window_coordinates)) {
            # If there are previous moving windows coordinates

            for (index in seq.int(
              from = ifelse(length(self$old_moving_window_coordinates) >= 3,
                length(self$old_moving_window_coordinates) - 2,
                1
              ),
              to = length(self$old_moving_window_coordinates))) {
              coords <- self$old_moving_window_coordinates[[index]]
              state_plot <- state_plot +
                annotate("rect",
                  xmin = coords[[1]][1],
                  xmax = coords[[2]][1],
                  ymin = coords[[1]][2],
                  ymax = coords[[2]][2],
                  color = "black",
                  alpha = .2
                )+
                annotate("label",
                  x = coords[[1]][1] + 0.25,
                  y = coords[[1]][2] + 0.25,
                  label = paste0(index)
                )
            }
          }

          if (!is.null(self$moving_window_coordinates)) {
            # We add the moving window on top of the plot
            coords <- self$moving_window_coordinates
            # coords[[1]] <- coords[[1]] - 0.25
            # coords[[2]] <- coords[[2]] + 0.25
            state_plot <- state_plot +
              annotate("rect",
                xmin = max(coords[[1]][1],1),
                xmax = min(coords[[2]][1],self$Q1_max),
                ymin = max(coords[[1]][2],1),
                ymax = min(coords[[2]][2],self$Q2_max),
                color = "red",
                alpha = .2
              )
          }


          state_plot <- state_plot +
            ggnewscale::new_scale_color() +
            geom_point(aes(
              x = Q1,
              y = Q2,
              size = BICL,
              color = isMaxBICL,
              alpha = BICL,
            )) +
          guides(color = guide_legend(title = "Is max value\nof BICL ?")) +
            ggnewscale::new_scale_color() +
            scale_colour_hue(l = 45, drop = FALSE) +
            geom_point(aes(
              x = Q1,
              y = Q2,
              color = clusteringComplete
            )) +
            guides(color = guide_legend(title = "Is the clustering complete ?"))
          ggtitle("State space for ", toString(self$net_id))


      print(state_plot)
      }
    },

    #' Function to greedily explore state of space looking for the mode
    #' and storing the models discovered along the way
    #'
    #' @param starting_point A vector of the two coordinates
    #' c(Q1,Q2) which are the starting point
    #' @param max_step_without_improvement defaults to 3, the
    #' number of steps to try improving before stopping the search
    #' @export
    #' @return c(Q1_mode, Q2_mode) which indicates the Q1 and Q2
    #' for which the BICL was maximal
    greedy_exploration = function(starting_point,
    max_step_without_improvement = 3) {
      # Initialize
      current_Q1 <- starting_point[1]
      current_Q2 <- starting_point[2]

      # Setting the starting point in the exploration order list
      self$exploration_order_list <- append(
        self$exploration_order_list,
        list(list(c(current_Q1, current_Q2)))
      )

      # Retrieve the index in the exploration order list
      index_in_exploration_order_list <- length(self$exploration_order_list)

      max_BICL_value <- -Inf
      max_BICL_coordinates <- NULL
      step_without_improvement <- 0
      max_BICL_has_improved <- TRUE
      step <- 0

      while (
        max_BICL_has_improved &&
        step < self$global_opts$Q1_max * self$global_opts$Q2_max
        ) {
        # The loop explores the space greedily
        if (self$global_opts$verbosity >= 4) {
          cat("\n----")
          cat(
            "\nExploring around Q = (",
            toString(c(current_Q1, current_Q2)), ")"
          )
        }

        # Appending the current point to the exploration order list
        self$exploration_order_list[[index_in_exploration_order_list]] <- append(
          self$exploration_order_list[[index_in_exploration_order_list]],
          list(c(current_Q1, current_Q2))
        )


        # The current model considered
        current_model <- self$model_list[[current_Q1, current_Q2]]

        neighbors <- list(c(1, 0), c(0, 1))
        # We loop through the neighbors of the current point
        for (neighbor in neighbors) {
          next_Q1 <- neighbor[1] + current_Q1
          next_Q2 <- neighbor[2] + current_Q2

          # Initialize
          next_Z_init <- vector("list", self$M)

          if (!self$point_is_in_limits(c(next_Q1, next_Q2))) {
            # The value is out of the allowed values, we quit this iteration
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\n(", toString(c(next_Q1, next_Q2)), ") is an INVALID neighbor for (",
                toString(c(current_Q1, current_Q2)), ").\nSkipping to next."
              )
            }
            next
          }

          # Otherwise the value is in the domain we want to explore
          if (self$global_opts$verbosity >= 4) {
            cat(
              "\n(", toString(c(next_Q1, next_Q2)),
              ") is a VALID neighbor for (",
              toString(c(current_Q1, current_Q2)), ")."
            )
          }

          # If we are splitting on the rows
          if (neighbor[[1]] == 1) {
            next_model <- self$split_clustering(current_model)
          }

          # If we are splitting on the columns
          if (neighbor[[2]] == 1) {
            next_model <- self$split_clustering(
              current_model,
              axis = "col"
            )
          }

          if (is.null(self$model_list[[next_Q1, next_Q2]])) {
            # If the point hasn't been seen yet the value is set
            self$model_list[[next_Q1, next_Q2]] <- next_model
          } else {
            # If the point has been seen
            # We append the min BICL to the discarded_model_list
            self$discarded_model_list[[next_Q1, next_Q2]] <- append(
              self$discarded_model_list[[next_Q1, next_Q2]], c(
                self$model_list[[next_Q1, next_Q2]],
                next_model
              )[[which.min(
                c(
                  self$model_list[[next_Q1, next_Q2]]$BICL,
                  next_model$BICL
                )
              )]]
            )

            # Performs the comparison and chooses the model that maximizes BICL
            self$model_list[[next_Q1, next_Q2]] <- c(
              self$model_list[[next_Q1, next_Q2]],
              next_model
            )[[which.max(
              c(
                self$model_list[[next_Q1, next_Q2]]$BICL,
                next_model$BICL
              )
            )]]
          }
        }

        # Now that all possible neighbors have been visited
        # we set the next point from which we'll loop again
        # by finding the neighbor with max BICL
        best_neighbor <- neighbors[[which.max(c(
          ifelse(current_Q1 + neighbors[[1]][[1]] <= self$global_opts$Q1_max && current_Q2 + neighbors[[1]][[2]] <= self$global_opts$Q2_max,
          self$model_list[[current_Q1 + neighbors[[1]][[1]],current_Q2 + neighbors[[1]][[2]]]]$BICL, # If the Q1 and Q2 are inbound we test them
          -Inf), # Else we compare to -Inf

          ifelse(current_Q1 + neighbors[[2]][[1]] <= self$global_opts$Q1_max && current_Q2 + neighbors[[2]][[2]] <= self$global_opts$Q2_max,
          self$model_list[[current_Q1 + neighbors[[2]][[1]], current_Q2 + neighbors[[2]][[2]]]]$BICL,
          -Inf
          )
        ))]]
        best_neighbor <- c(current_Q1 + best_neighbor[1], current_Q2 + best_neighbor[2])

        # Now we set the current Q to the best neighbor
        current_Q1 <- best_neighbor[1]
        current_Q2 <- best_neighbor[2]

        if (self$model_list[[best_neighbor[1], best_neighbor[2]]]$BICL > max_BICL_value){
          # If the neighbor we found is best than
          # the previous mode we update and go for one more iteration
          max_BICL_value <- self$model_list[[best_neighbor[1], best_neighbor[2]]]$BICL
          max_BICL_coordinates <- c(
            best_neighbor[1],
            best_neighbor[2]
          )
          step_without_improvement <- 0
          max_BICL_has_improved <- TRUE
        } else {
          # Else we've found a local mode
          step_without_improvement <- step_without_improvement + 1
          if (step_without_improvement >= max_step_without_improvement) {
            max_BICL_has_improved <- FALSE
          }
        }

        if (max_BICL_has_improved && step_without_improvement == 0) {
          # If the BICL improved in this round
          end_of_text <- " and the BICL improved at this step. Going for another step"
        } else if (max_BICL_has_improved && step_without_improvement <= max_step_without_improvement) {
          end_of_text <- paste0(" and the BICL hasn't improved at this step. ", step_without_improvement, "/", max_step_without_improvement, " steps without improvement.")
        } else {
          end_of_text <- paste0(" and the BICL hasn't improved ", max_step_without_improvement, " times in a row, stopping search.")
        }

        if (self$global_opts$verbosity >= 4) {
          cat("\nFor this round the best neighbor is: ", toString(best_neighbor), end_of_text)
        }

        # We increase the step
        step <- step + 1
      }

      self$exploration_order_list[[index_in_exploration_order_list]] <- self$exploration_order_list[[index_in_exploration_order_list]][-1]

      # Return the coordinates of the max BICL that has been found
      return(max_BICL_coordinates)
    },

    #' Optimization from a given Z_init
    #' 
    #' The functions takes no parameters
    #' 
    #' @return nothing but stores the models in the model list
    optimize_from_zinit = function() {
      model_list <- vector(
        "list",
        self$global_opts$Q1_max * self$global_opts$Q2_max
      )
      dim(model_list) <- c(self$global_opts$Q1_max, self$global_opts$Q2_max)

      #' The function fit a fitBipartiteSBMPop object with the wanted
      #' parameters
      #'
      #' @param q1 the row cluster number
      #' @param q2 the col cluster number
      #' @param Z the given clustering to fit from
      #' @param Cpi (optional) the supports if provided
      #' @param Calpha (optional) the alpha support, if provided
      #'
      #' @return a fitBipartite fitted object
      optimize_init <- function(q1, q2, Z, Cpi = NULL, Calpha= NULL) {
        fit <- fitBipartiteSBMPop$new(
          A = self$A,
          Q = c(q1, q2),
          free_mixture_row = self$free_mixture_row,
          free_mixture_col = self$free_mixture_col,
          free_density = self$free_density,
          Cpi = Cpi,
          Calpha = Calpha,
          init_method = "given",
          distribution = self$distribution,
          Z = Z, # This contains the Z_init to fit
          net_id = self$net_id,
          fit_opts = self$fit_opts,
        )
        if (self$global_opts$verbosity >=4) {
          cat("\nOptimizing from Z_init for", toString(c(q1,q2)))
        }

        fit$optimize()
        return(fit)
      }

      # We fill the model_list from the Z_init provided
      # assuming that Z_init is a bi-dimensional list
      for (q1 in seq.int(self$global_opts$Q1_max)) {
        for (q2 in seq.int(self$global_opts$Q2_max)) {
          if (!is.null(self$Z_init[[q1, q2]])) {
            model_list[[q1, q2]] <- optimize_init(
              q1, q2, self$Z_init[[q1, q2]]
              # FIXME we dont provide Cpi nor Calpha, thus the fitted points
              # will have iid-like fits.
              # This will be overwritten when the moving window will perform
            )
          }
        }
      }
      # FIXME This doesn't fill the model_list
      # lapply(
      #   seq.int(self$global_opts$Q1_max),
      #   function(q1) {
      #     lapply(
      #       seq.int(self$global_opts$Q2_max),
      #       function(q2) {
      #         if (!is.null(self$Z_init[[q1, q2]])) {
      #           model_list[[q1, q2]] <- optimize_init(
      #             q1, q2, self$Z_init[[q1, q2]]
      #             # FIXME we dont provide Cpi nor Calpha, thus the fitted points
      #             # will have iid-like fits.
      #             # This will be overwritten when the moving window will perform
      #           )
      #         }
      #       }
      #     )
      #   }
      # )
      # model_list contains the fitted objects and we can set the self$
      self$model_list <- model_list
      dim(model_list) <- c(self$global_opts$Q1_max, self$global_opts$Q2_max)
      self$store_criteria_and_best_fit()
    },

    #' Burn-in method to start exploring the state of space
    #'
    #' The functions takes no parameters but modify the object
    #'
    #' @return nothing; but stores the values
    burn_in = function() {
      start_time <- Sys.time()
      # The function fit M fitBipartite from spectral clust for Q = (1,2) and Q = (2,1)
      # DONE : implement the M fitBipartite from spectral clust


      if (self$global_opts$verbosity >= 2) {
        cat("\n==== Beginning Burn in ====")
      }
      if (!is.null(self$Z_init)) {
        # If a Z_init was provided we perform its optimization to save time
        if (self$global_opts$verbosity >= 2) {
          cat("\nA Z_init was provided, fitting the provided points.")
        }
        self$optimize_from_zinit()
      } else { # TODO I can choose to fit from Z_init AND THEN perform greedy
               # exploration
        # With no Z_init we perform a greedy exploration
        self$separated_inits <- vector("list", 4) # The first coordinate is Q1, the second is Q2
        dim(self$separated_inits) <- c(2,2)
        if (self$global_opts$verbosity >= 2) {
          cat("\nFitting good initilization points.")
        }

        # Init for Q = (1,2)
        if (self$global_opts$verbosity >= 4) {
          cat("\nFitting ", self$M, " networks for Q = (", toString(c(1, 2)), ")\n")
        }

        self$separated_inits[[1,2]] <- lapply(
          seq.int(self$M),
          function(m) {
            fitBipartiteSBMPop$new(
              A = list(self$A[[m]]),
              Q = c(1, 2),
              free_mixture_row = FALSE, # There can't be free mixture with 1 net
              free_mixture_col = FALSE, # There can't be free mixture with 1 net
              free_density = self$free_density,
              distribution = self$distribution,
              init_method = "spectral",
              fit_opts = self$fit_opts
            )
          }
            )

        # Fitting each model
        lapply(
          seq.int(self$M),
          function(m) {
            self$separated_inits[[1,2]][[m]]$optimize()
          }
        )

        if (self$global_opts$verbosity >= 4) {
          sep_vbounds <- lapply(seq.int(self$M),
          function(m){
            max(self$separated_inits[[1,2]][[m]]$vbound)
          })
          cat(
            "Finished fitting ", self$M, " networks for Q = (", toString(c(1, 2)),
            ")\nSeparated Variational Bounds for the networks", toString(seq.int(self$M)),
            ":\n", toString(sep_vbounds), "\n"
          )
        }

        # Init for Q = (2,1)
        if (self$global_opts$verbosity >= 4) {
          cat("\nFitting ", self$M, " networks for Q = (", toString(c(2, 1)), ")\n")
        }

        self$separated_inits[[2,1]] <- lapply(
          seq.int(self$M),
          function(m) {
            fitBipartiteSBMPop$new(
              A = list(self$A[[m]]), 
              Q = c(2, 1),
              free_mixture_row = FALSE, # There can't be free mixture with 1 net
              free_mixture_col = FALSE, # There can't be free mixture with 1 net
              free_density = self$free_density,
              distribution = self$distribution,
              init_method = "spectral",
              fit_opts = self$fit_opts
            )
          }
        )

        # Fitting each model
        lapply(
          seq.int(self$M),
          function(m) {
            self$separated_inits[[2,1]][[m]]$optimize()
          }
        )
      if (self$global_opts$verbosity >= 4) {
        sep_vbounds <- lapply(
          seq.int(self$M),
          function(m) {
            max(self$separated_inits[[2,1]][[m]]$vbound)
          }
        )
        cat(
          "Finished fitting ", self$M, " networks for Q = (", toString(c(2, 1)),
          ")\nSeparated Variational Bounds for the networks", toString(seq.int(self$M)),
          ":\n", toString(sep_vbounds), "\n"
        )
      }

        # Here we match the clusters from the M fit objects
        # By using the order of the marginal laws
        if (self$global_opts$verbosity >= 4) {
          cat("\nBeginning to match results for the Separated BiSBMs.")
        }

        for (m in seq.int(self$M)) {
          current_m_init <- self$separated_inits[[1,2]][[m]]

          # The clustering are one hot encoded because the permutations are
          # easier to perform

          # One hot encoded row (1 cluster)
          row_clustering <- .one_hot(current_m_init$Z[[1]][[1]], 1)
          # One hot encoded cols (2 clusters)
          col_clustering <- .one_hot(current_m_init$Z[[1]][[2]], 2)

          # The row clustering are reordered according to their marginal distribution
          prob1 <- as.vector(
            current_m_init$pi[[1]][[2]] %*% t(current_m_init$MAP$alpha)
          )

          p1 <- order(prob1)

          row_clustering <- row_clustering[, p1]

          # The col clustering are reordered according to their marginal distribution
          prob2 <- as.vector(
            current_m_init$pi[[1]][[1]] %*% current_m_init$MAP$alpha
          )
          p2 <- order(prob2)
          col_clustering <- col_clustering[, p2]


          # The clustering are reverse one hot encoded
          row_clustering <- .rev_one_hot(row_clustering)
          col_clustering <- .rev_one_hot(col_clustering)

          self$separated_inits[[1,2]][[m]]$Z[[1]][[1]] <- row_clustering
          self$separated_inits[[1,2]][[m]]$Z[[1]][[2]] <- col_clustering

          if (self$global_opts$verbosity >= 4) {
            cat(
              "\nNetwork:", m,
              "\nOrder lines:", toString(p1),
              "\nOrder columns:", toString(p2),
              "\n"
            )
          }
        }

        if (self$global_opts$verbosity >= 4) {
          cat("\nMatching finished.\n")
          cat("Beginning to combine networks.")
        }

        if (self$global_opts$verbosity >= 4) {
          cat("\nFitting full model Q=(1,1) for ", self$M)
        }

        self$model_list[[1, 1]] <- fitBipartiteSBMPop$new(
          A = self$A,
          Q = c(1, 1),
          free_mixture_row = FALSE, # There can't be free mixture with 1 cluster
          free_mixture_col = FALSE, # There can't be free mixture with 1 cluster
          free_density = self$free_density,
          fit_opts = self$fit_opts,
          distribution = self$distribution,
          greedy_exploration_starting_point = c(1,1),
        )

        self$model_list[[1,1]]$optimize()

        # Here we combine the networks to fit a
        # fitBipartite object on the M networks

        # We retrieve the clustering for the M (1,2) separated models
        M_clusterings_1_2 <- lapply(
          seq.int(self$M),
          function(m) {
            # We add [[1]] after Z because we fitted
            # only one network with the objects stored
            # where the class is supposed to store more
            self$separated_inits[[1,2]][[m]]$Z[[1]]
          }
        )

        # We retrieve the clustering for the M (2,1) separated models
        M_clusterings_2_1 <- lapply(
          seq.int(self$M),
          function(m) {
            # We add [[1]] after Z because we fitted
            # only one network with the objects stored
            # where the class is supposed to store more
            self$separated_inits[[2,1]][[m]]$Z[[1]]
          }
        )

        self$separated_inits[[1,2]] <- fitBipartiteSBMPop$new(
          A = self$A, Q = c(1, 2),
          free_mixture_row = FALSE, # There can't be free mixture with 1 cluster
          free_mixture_col = self$free_mixture_col,
          free_density = self$free_density,
          Z = M_clusterings_1_2,
          init_method = "given",
          distribution = self$distribution,
          fit_opts = self$fit_opts,
          greedy_exploration_starting_point = c(1,2)
        )

        self$separated_inits[[2,1]] <- fitBipartiteSBMPop$new(
          A = self$A, Q = c(2, 1),
          free_mixture_row = self$free_mixture_row,
          free_mixture_col = FALSE, # There can't be free mixture with 1 cluster
          free_density = self$free_density,
          Z = M_clusterings_2_1,
          init_method = "given",
          distribution = self$distribution,
          fit_opts = self$fit_opts,
          greedy_exploration_starting_point = c(2,1)
        )

        if (self$global_opts$verbosity >= 4) {
          cat("\nFitting the combined colBiSBMs.")
        }
        # Here we fit the models
        lapply(seq.int(2),
        function(index){
          # The index used here are a little trick to allow the use
          if (self$global_opts$verbosity >= 4){
            cat(
              "\nFitting the ", self$M, " networks for Q = (",
              toString(c(index, 3 - index)), ")."
            )
          }

          self$separated_inits[[index,3 - index]]$optimize()
        })

        if (self$global_opts$verbosity >= 4) {
          cat("\nFinished fitting the colBiSBM.")
          cat("\nResults for the the points :")
          for (index in c(1,2)){
            cat("\nQ = (", toString(c(index, 3 - index)), ") :")
            cat(
              "\n\tvbound:",
              toString(self$separated_inits[[index,3 - index]]$vbound)
            )
            cat(
              "\n\tICL:",
              toString(self$separated_inits[[index,3 - index]]$ICL)
            )
            cat(
              "\n\tBICL:",
              toString(self$separated_inits[[index,3 - index]]$BICL)
            )
          }
        }

        # Store the given initialization
        self$model_list[1, 2] <- self$separated_inits[1, 2]
        self$model_list[2, 1] <- self$separated_inits[2, 1]

        # We go looking for the mode with a greedy approach
        # Visiting each of the neighbors
        if (self$global_opts$verbosity >= 3) {
          cat("\n Greedy exploration to find a first mode.")
        }

        # Greedy exploration from (1,2)
        mode_1_2 <- self$greedy_exploration(c(1, 2))
        if (self$global_opts$verbosity >= 4){
          cat("\nFrom (", toString(c(1,2)),") the mode is at: (", toString(mode_1_2),").")
        }

        # Greedy exploration from (2,1)
        mode_2_1 <- self$greedy_exploration(c(2, 1))
        if (self$global_opts$verbosity >= 4) {
          cat("\nFrom (", toString(c(2, 1)), ") the mode is at: (", toString(mode_2_1), ").")
          # Plot the state of space and it's exploration
          self$state_space_plot()
        }

        # Finding the best of the two modes
        best_mode <- list(mode_1_2, mode_2_1)[[which.max(c(
          self$model_list[[mode_1_2[1], mode_1_2[2]]]$BICL,
          self$model_list[[mode_2_1[1], mode_2_1[2]]]$BICL
        ))]]
      }

      self$store_criteria_and_best_fit()

      if(self$global_opts$verbosity >= 2) {
        cat(
          "\n==== Finished Burn in",
          " for networks ", self$net_id, " in ",
          format(Sys.time() - start_time, digits = 3),
          " ===="
        )
      }

      if(self$global_opts$verbosity >= 2) {
      self$print_metrics()

      }
    },

    optimize = function() {
      # The burn_in step computes models with a greedy approach
      self$burn_in()
      improved <- TRUE
      nb_pass <- 0
      tolerance <- 10e-3
      self$store_criteria_and_best_fit()
      Q <- which(self$BICL == max(self$BICL), arr.ind = TRUE)
      current_max_BICL <- self$model_list[[Q[1], Q[2]]]$BICL

      if (self$global_opts$verbosity >= 2) {
        cat("\n==== Beginning the passes of moving windows ====")
      }

      self$global_opts$nb_models <- ceiling(self$global_opts$nb_models/2)
      while (improved & nb_pass < self$global_opts$max_pass) {
        if (self$global_opts$verbosity >= 2) {
          cat(
            "\n==== Starting pass number ", nb_pass + 1,
            " for networks ", self$net_id, " ====\n"
          )
        }

        # Storing the current important values
        current_max_BICL <- self$model_list[[Q[1],Q[2]]]$BICL

        # Perform another iteration of the moving window
        self$moving_window(Q, depth = self$global_opts$depth)

        # Improvement criterion
        has_the_mode_moved <- all(Q != self$best_fit$Q)
        self$improved <- (self$best_fit$BICL - current_max_BICL) > tolerance ||
          has_the_mode_moved

        # We now set the new Q and check if the fit is better
        Q <- which(self$BICL == max(self$BICL), arr.ind = TRUE)
        improved <- self$improved
        nb_pass <- nb_pass + 1
        self$truncate_discarded_model_list()

        if (self$global_opts$verbosity >= 2) {
          self$print_metrics()
        }
      }

      if (self$global_opts$verbosity >= 2) {
        cat(
          "\n==== Finished the passes of moving windows ====",
          "\n==== Computing separated BiSBM BICL ===="
        )
      }

      # Now compare to the sepBiSBM
      if (self$M > 1) {
        # If there is more than one network we launch the choose procedure
        self$choose_joint_or_separated()
      } else {
        # Otherwise we're inside a choose procedure and store the BICL
        self$sep_BiSBM$BICL <- c(self$best_fit$BICL)
      }

    },

    #' The moving window application
    #'
    #' @description
    #' This method is a moving windows
    #' over the Q1xQ2 space for the number of clusters
    #' @noMd
    #' @noRd
    #' @param center is the coordinates (Q1, Q2) in the model list of the mode
    #' @param depth is how far away from the center
    #' the function should be applied in a grid style
    #' going from center - (depth,depth) to center + (depth, depth)
    #' @return nothing; but updates the object by adding new models
    moving_window = function(center, depth = 1) {
      # Each split & merge can be parallelized from points
      # But need to be finished when comparing the BICL

      Q1_mode <- center[[1]]
      Q2_mode <- center[[2]]

      # Checking if the window's bound is in domain
      if (!self$point_is_in_limits(c(Q1_mode - depth, Q2_mode - depth)) ||
      !self$point_is_in_limits(c(Q1_mode + depth, Q2_mode + depth))) {
        warning(paste0(
          "\nThe windows is (partially) out of domain !",
          "\nTrying to go from (", toString(c(Q1_mode - depth, Q2_mode - depth)),
          ") to (", toString(c(Q1_mode + depth, Q2_mode + depth)), ").",
          "\nMax windows should be :",
          "\n(", toString(c(1, self$global_opts$Q2_max)),
          ")---(", toString(c(self$global_opts$Q1_max, self$global_opts$Q2_max)),
          ")",
          "\n  ||       || ",
          "\n(", toString(c(1, 1)),
          ")---(",
          toString(c(self$global_opts$Q1_max, 1)),
          ")",
          "\nThe window will only work on valid configurations"
        ))

      }

      if (self$global_opts$verbosity >= 4) {
        cat(
          "\nMoving window around (",
          toString(center), ") with a depth of ", depth
        )
      }

      if (!is.null(self$moving_window_coordinates)) {
        # If there was already a previous iteration of the moving window
        if (any(
          self$moving_window_coordinates[[1]] != c(Q1_mode - depth, Q2_mode - depth)
        ) &&
          any(
            self$moving_window_coordinates[[2]] != c(Q1_mode + depth, Q2_mode + depth)
          )) {
          # If the moving window is not on the same point
          # its coords are appended as old coords
          self$old_moving_window_coordinates <- append(
            self$old_moving_window_coordinates,
            list(self$moving_window_coordinates)
          )
        }
      }

      self$moving_window_coordinates <- list(
        c(Q1_mode - depth, Q2_mode - depth),
        c(Q1_mode + depth, Q2_mode + depth)
      )

      # We perform the Forward Pass
      if (self$global_opts$verbosity >= 4) {
        cat("\nBeginning the forward pass.")
      }
      # Forward pass, where we split

      for (current_Q1 in seq.int(from = Q1_mode -depth, to = Q1_mode + depth)) {
        if (current_Q1 < 1 || current_Q1 > self$global_opts$Q1_max) {
          # Discard the value if it's out of bounds
          next
        }
        for (current_Q2 in seq.int(from = Q2_mode-depth, to = Q2_mode + depth)){
          if (current_Q2 < 1 || current_Q2 > self$global_opts$Q2_max) {
            # Discard the value if it's out of bounds
            next
          }

          current_model_Q <- c(current_Q1, current_Q2)

          # We list the split origins for the model and we'll keep the best
          # BICL and store the other one in the discarded model list
          wanted_model_different_splits_origin <- list()
          left_model <- NULL
          bottom_model <- NULL

          if (self$global_opts$verbosity >= 4) {
            cat(
              "\n\nLooking for possible origins by splitting to Q = (",
              toString(current_model_Q), ")."
            )
          }

          # If the wanted model already exists in the model_list we store it as a possible model
          if (self$point_is_in_limits(current_model_Q) &&
          !is.null(self$model_list[[current_model_Q[1], current_model_Q[2]]])) {
            if (self$global_opts$verbosity >= 4) {
              cat("\nA model was already fitted here ! Storing it to compare")
            }

            wanted_model_different_splits_origin <-
              append(
                wanted_model_different_splits_origin,
                self$model_list[[current_model_Q[1], current_model_Q[2]]]
              )
          }else if (self$point_is_in_limits(current_model_Q) &&
          is.null(self$model_list[[current_model_Q[1], current_model_Q[2]]]) &&
          self$global_opts$verbosity >= 4) {
            cat("\nNo model already fitted for the point")
          }

          left_model_Q <- current_model_Q + c(-1,0)

          # Checking if the left neighbor exists
          if (self$point_is_in_limits(left_model_Q) &&
            !is.null(self$model_list[[left_model_Q[1], left_model_Q[2]]])) {
            # If the left neighbor exist then we can split from it
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nThe left neighbor of (",
                toString(current_model_Q),
                ") exists. It is (",
                toString(left_model_Q),
                ").\nFitting the possible row splits from it."
              )
            }

            left_model <- self$model_list[[left_model_Q[1], left_model_Q[2]]]

            wanted_model_different_splits_origin <- append(
              wanted_model_different_splits_origin,
              self$split_clustering(left_model)
            )
          } else if (self$point_is_in_limits(left_model_Q) &&
            is.null(self$model_list[[left_model_Q[1], left_model_Q[2]]]) &&
            self$global_opts$verbosity >= 4) {
              cat("\nNo left neighbor already fitted")
            }

          bottom_model_Q <- current_model_Q + c(0,-1)

          # Checking if the bottom neighbor exists
          if (self$point_is_in_limits(bottom_model_Q)&&
          !is.null(self$model_list[[bottom_model_Q[1], bottom_model_Q[2]]])) {
            # If the bottom neighbor exist then we can split from it
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nThe bottom neighbor of (",
                toString(current_model_Q),
                ") exists. It is (",
                toString(
                  bottom_model_Q),
                  ").\nFitting the possible column splits from it."
              )
            }

            bottom_model <- self$model_list[[bottom_model_Q[1], bottom_model_Q[2]]]

            wanted_model_different_splits_origin <- append(
              wanted_model_different_splits_origin,
              self$split_clustering(bottom_model, axis = "col")
            )
          } else if (self$point_is_in_limits(bottom_model_Q) &&
          is.null(self$model_list[[bottom_model_Q[1], bottom_model_Q[2]]]) &&
          self$global_opts$verbosity >= 4) {
            cat("\nNo bottom neighbor already fitted")
          }

          # If the point has no predecessor or current model
          # a spectral is fitted
          if (length(wanted_model_different_splits_origin) == 0) {
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nNo possible origins nor fitted model for the point Q=(",
                toString(current_model_Q),
                "). Fitting a spectral at this point."
              )
            }
            # OPTIONAL TODO : release the if, to test with more
            # inits in the grid
            spectral_init <- fitBipartiteSBMPop$new(
              A = self$A,
              Q = current_model_Q,
              free_mixture_row = self$free_mixture_row,
              free_mixture_col = self$free_mixture_col,
              free_density = self$free_density,
              init_method = "spectral",
              distribution = self$distribution,
              net_id = self$net_id,
              fit_opts = self$fit_opts
            )
            spectral_init$optimize()
            wanted_model_different_splits_origin <- append(
              wanted_model_different_splits_origin,
              spectral_init
            )
          }

          # Now we have the different models from different origins
          # and we can select the one that maximizes the BICL

          # We compute an intermediate list of the BICLs
          wanted_model_different_splits_origin_BICL <- NULL
          wanted_model_different_splits_origin_BICL <- lapply(
            seq_along(wanted_model_different_splits_origin),
            function(origin) {
              wanted_model_different_splits_origin[[origin]]$BICL
            }
          )
          # Adding the best of the possible models to the model_list
          self$model_list[[current_model_Q[1], current_model_Q[2]]] <-
            wanted_model_different_splits_origin[[which.max(wanted_model_different_splits_origin_BICL)]]

          if (self$global_opts$verbosity >= 4) {
            cat(
              "\nThe preferred origin for (",
              toString(current_model_Q),
              ") is the Q = (",
              toString(wanted_model_different_splits_origin[[which.max(wanted_model_different_splits_origin_BICL)]]$Q),
              ") model."
            )
          }

          # Adding the other possible models to the discarded_model_list
          self$discarded_model_list[[current_model_Q[1], current_model_Q[2]]] <-
          append(
            self$discarded_model_list[[current_model_Q[1], current_model_Q[2]]],
            wanted_model_different_splits_origin[-which.max(wanted_model_different_splits_origin_BICL)]
          )
        }
        }

        self$state_space_plot()

      if (self$global_opts$verbosity >= 4) {
        cat("\nEnd of the Forward pass.\n")
      }

      # Backward pass, where we merge
      for (current_Q1 in seq.int(from = Q1_mode + depth, to = Q1_mode - depth)) {
        if (current_Q1 < 1 || current_Q1 > self$global_opts$Q1_max) {
          next
        }
        for (current_Q2 in seq.int(from = Q2_mode + depth, to = Q2_mode - depth)) {
          if (current_Q2 < 1 || current_Q2 > self$global_opts$Q2_max) {
            next
          }
          # Current model to merge
          current_model_Q <- c(current_Q1, current_Q2)

          # We list the split origins for the model and we'll keep the best
          # BICL and store the other one in the discarded model list
          wanted_model_different_merges_origin <- list()
          right_model <- NULL
          top_model <- NULL

          if (self$global_opts$verbosity >= 4) {
            cat(
              "\n\nLooking for possible origins by merging to Q = (",
              toString(current_model_Q), ")."
            )
          }

          # If the wanted model already exists in the model_list
          # we store it as a possible model
          if (self$point_is_in_limits(current_model_Q) &&
            !is.null(self$model_list[[current_model_Q[1], current_model_Q[2]]])) {
            if (self$global_opts$verbosity >= 4) {
              cat("\nA model was already fitted here ! Storing it to compare")
            }

            wanted_model_different_merges_origin <-
              append(
                wanted_model_different_merges_origin,
                self$model_list[[current_model_Q[1],current_model_Q[2]]]
              )
          } else if (self$point_is_in_limits(current_model_Q) &&
            is.null(self$model_list[[current_model_Q[1], current_model_Q[2]]])&&
            self$global_opts$verbosity >= 4) {
            cat("\nNo model already fitted for the point")
          }

          # Checking if the right neighbor exists

          right_model_Q <- current_model_Q + c(1,0)

          if (self$point_is_in_limits(right_model_Q) &&
            !is.null(self$model_list[[right_model_Q[1], right_model_Q[2]]])) {
            # If the right neighbor exist then we can split from it
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nThe right neighbor of (",
                toString(current_model_Q),
                ") exists. It is (",
                toString(right_model_Q),
                ").\nFitting the possible row merges from it."
              )
            }

            right_model <- self$model_list[[right_model_Q[1], right_model_Q[2]]]

            wanted_model_different_merges_origin <- append(
              wanted_model_different_merges_origin,
              self$merge_clustering(right_model)
            )
          } else if (self$point_is_in_limits(right_model_Q) &&
            is.null(self$model_list[[right_model_Q[1], right_model_Q[2]]]) &&
            self$global_opts$verbosity >= 4) {
            cat("\nNo right neighbor already fitted")
          }

          # Checking if the top neighbor exists
          top_model_Q <- current_model_Q + c(0, 1)

          if (self$point_is_in_limits(top_model_Q) &&
            !is.null(self$model_list[[top_model_Q[1], top_model_Q[2]]])) {
            # If the top neighbor exist then we can split from it
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nThe top neighbor of (",
                toString(current_model_Q),
                ") exists. It is (",
                toString(
                  top_model_Q
                ),
                ").\nFitting the possible column merges from it."
              )
            }

            top_model <- self$model_list[[top_model_Q[1],top_model_Q[2]]]

            wanted_model_different_merges_origin <- append(
              wanted_model_different_merges_origin,
              self$merge_clustering(top_model, axis = "col")
            )
          } else if (self$point_is_in_limits(top_model_Q) &&
            is.null(self$model_list[[top_model_Q[1], top_model_Q[2]]]) &&
            self$global_opts$verbosity >= 4) {
            cat("\nNo top neighbor already fitted")
          }

          # If the point has no predecessor or current model
          # a spectral is fitted
          if (length(wanted_model_different_merges_origin) == 0) {
            if (self$global_opts$verbosity >= 4) {
              cat(
                "\nNo possible origins nor fitted model for the point Q=(",
                toString(current_model_Q),
                "). Fitting a spectral at this point."
              )
            }
            # OPTIONAL TODO : release the if, to test with more
            # inits in the grid
            spectral_init <- fitBipartiteSBMPop$new(
              A = self$A,
              Q = current_model_Q,
              free_mixture_row = self$free_mixture_row,
              free_mixture_col = self$free_mixture_col,
              free_density = self$free_density,
              init_method = "spectral",
              distribution = self$distribution,
              net_id = self$net_id,
              fit_opts = self$fit_opts
            )
            spectral_init$optimize()
            wanted_model_different_merges_origin <- append(
              wanted_model_different_merges_origin,
              spectral_init
            )
          }

          # Now we have the different models from different origins
          # and we can select the one that maximizes the BICL

          # We compute an intermediate list of the BICLs
          wanted_model_different_merges_origin_BICL <- NULL
          wanted_model_different_merges_origin_BICL <- lapply(
            seq_along(wanted_model_different_merges_origin),
            function(origin) {
              wanted_model_different_merges_origin[[origin]]$BICL
            }
          )
          # Adding the best of the possible models to the model_list
          self$model_list[[current_model_Q[1], current_model_Q[2]]] <-
            wanted_model_different_merges_origin[[which.max(wanted_model_different_merges_origin_BICL)]]

          if (self$global_opts$verbosity >= 4) {
            cat(
              "\nThe preferred origin for (",
              toString(current_model_Q),
              ") is the Q = (",
              toString(
                wanted_model_different_merges_origin[[which.max(wanted_model_different_merges_origin_BICL)]]$Q
                ),
              ") model."
            )
          }

          # Adding the other possible models to the discarded_model_list
          self$discarded_model_list[[current_model_Q[1], current_model_Q[2]]] <- append(
            self$discarded_model_list[[current_model_Q[1], current_model_Q[2]]],
            wanted_model_different_merges_origin[-which.max(wanted_model_different_merges_origin_BICL)]
          )
        }
      }

      self$state_space_plot()

      if (self$global_opts$verbosity >= 4) {
        cat("\nEnd of the Backward pass.\n")
      }

      # After the two passes
      self$store_criteria_and_best_fit()

      if (self$global_opts$verbosity >= 4) {
        cat(
          "\nAfter the moving window around (",
          toString(c(Q1_mode, Q2_mode)),
          ") the mode is at : (",
          toString(self$best_fit$Q),
          ")."
        )
      }
    },

    truncate_discarded_model_list = function() {
      for (q1 in seq.int(self$global_opts$Q1_max)) {
        for (q2 in seq.int(self$global_opts$Q2_max)) {
          if (!is.null(self$discarded_model_list[[q1,q2]]) &&
          length(self$discarded_model_list[[q1,q2]]) > self$global_opts$nb_models - 1) {
            models_BICL <-
              sapply(seq_along(self$discarded_model_list[[q1, q2]]), function(s) {
                self$discarded_model_list[[q1, q2]][[s]]$BICL
              })
            self$discarded_model_list[[q1,q2]] <-
              self$discarded_model_list[[q1,q2]][order(models_BICL, decreasing = TRUE)][1:(self$global_opts$nb_models - 1)]
          }
        }
      }
    },

    point_is_in_limits = function(point) {
      Q1 <- point[[1]]
      Q2 <- point[[2]]

      return(Q1 > 0 && Q2 > 0 &&
        Q1 <= self$global_opts$Q1_max &&
        Q2 <= self$global_opts$Q2_max)
    },

    store_criteria_and_best_fit = function() {
      # Store the vbound, ICL and BICL into the appropriate lists

      lapply(seq.int(self$global_opts$Q1_max), function(q1) {
        lapply(seq.int(self$global_opts$Q2_max), function(q2) {
          current_model <- self$model_list[[q1, q2]]
          # If the model has been detected by the exploration
          if (!is.null(current_model)) {
            # The below expression handles the case where vbound is a null list
            self$vbound[q1, q2] <- ifelse(is.null(unlist(tail(self$model_list[[q1, q2]]$vbound, n = 1))),
              -Inf,
              unlist(tail(self$model_list[[q1, q2]]$vbound, n = 1))
            )
            self$ICL[q1, q2] <- current_model$ICL
            self$BICL[q1, q2] <- current_model$BICL
          }
        })
      })

      # Assign the best_fit
      self$best_fit <- self$model_list[[which.max(self$BICL)]]
    },

    compute_sep_BiSBM_BICL = function() {
      # Computes the sepBiSBM ICL to compare with the model
      # TODO See if I can parallelize
      self$sep_BiSBM$models <- lapply(seq.int(self$M), function(m) {
        sep_BiSBM <- bisbmpop$new(
          netlist = list(self$A[[m]]),
          distribution = self$distribution,
          free_mixture_row = FALSE,
          free_mixture_col = FALSE,
          free_density = FALSE,
          global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = self$global_opts$nb_cores
          )
        )
        sep_BiSBM$optimize()
        sep_BiSBM$best_fit
      })
      
      self$sep_BiSBM$BICL <- sapply(seq.int(self$M), function(m) {
        self$sep_BiSBM$models[[m]]$BICL  # We retrieve all the BICLs for
                                                  # separated models
      })

      self$sep_BiSBM$Z <- lapply(seq.int(self$M), function(m) {
        list(
          self$sep_BiSBM$models[[m]]$Z[[1]][[1]],
          self$sep_BiSBM$models[[m]]$Z[[1]][[2]]
        )
      })

      if (self$global_opts$verbosity >= 2) {
        cat("\n==== Finished fitting ", self$M, "sepBiSBM ====")
        cat("\n Total sepBiSBM BICL : ", sum(self$sep_BiSBM$BICL),"\n")
      }

    },

    choose_joint_or_separated = function() {
      # Now compare to the sepBiSBM
      if (is.null(self$sep_BiSBM$BICL)) {
        self$compute_sep_BiSBM_BICL()
      }

      if (self$global_opts$verbosity >= 1) {
        if (self$M > 1) {
          cat("\n==== Best fits criterion for the", self$M, "networks ====")
          cat("\nSep BiSBM total BICL: ", sum(self$sep_BiSBM$BICL))
          cat("\ncolBiSBM BICL:", self$best_fit$BICL)
          if (sum(self$sep_BiSBM$BICL) > self$best_fit$BICL) {
            cat("\nSeparated modelisation preferred.\n")
          } else {
            cat(
              "\nJoint modelisation preferred. With Q = (",
              toString(self$best_fit$Q),
              ").\n"
            )
          }
        } else {
          cat("\n==== Best fits criterion for the network ====")
          cat("\nOnly 1 network")
          cat("\nSep BiSBM total BICL: ", sum(self$sep_BiSBM$BICL))
        }
      }
    },

    print_metrics = function() {
      # Capturing the pretty print of matrices
      vbound_print <- paste0(capture.output(
        print(round(self$vbound))
      ), collapse = "\n")

      ICL_print <- paste0(capture.output(
        print(round(self$ICL))
      ), collapse = "\n")

      BICL_print <- paste0(capture.output(
        print(round(self$BICL))
      ), collapse = "\n")

      cat("\nvbound : \n", vbound_print, "\n\n")
      cat("ICL    : \n", ICL_print, "\n\n")
      cat("BICL   : \n", BICL_print, "\n\n")
      cat("Best fit at Q=(", toString(self$best_fit$Q), ")\n")
    }
  )
)

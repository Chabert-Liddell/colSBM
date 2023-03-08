#' An R6 Class object, a collection of model for population of LBM netowrks
#'
#' @export

lbmpop <- R6::R6Class(
  "lbmpop",
  #
  public = list(
    nr = NULL,
    nc = NULL,
    A = NULL,
    M = NULL,
    mask = NULL, # 1 for NA and 0 for observed
    model = NULL,
    net_id = NULL,
    model_list = NULL, # A list of size Q1max, of lists of size Q2max containing the models
    global_opts = NULL,
    fit_opts = NULL,
    fit_sbm = NULL,
    separated_inits = NULL, # A nested list : Q1 init containing Q2 init
                            # with each entry containing list of size M 
                            # storing the separated inits for the class
    Z_init = NULL,
    free_density = NULL,
    free_mixture = NULL,
    ICL_sbm = NULL,
    ICL = NULL,
    BICL = NULL,
    vbound = NULL,
    best_fit = NULL,
    logfactA = NULL,
    improved = list(forward = TRUE,
                    backward = TRUE),



    #' @description
    #' Create a new instance of the lbmpop object
    #' 
    #' This class is generally called via the user function # FIXME put the user function name
    #'
    initialize = function(netlist = NULL,
                          net_id = NULL,
                          model = "bernoulli",
                          free_density = FALSE,
                          free_mixture = FALSE,
                          fit_sbm = NULL,
                          Z_init = NULL,
                          global_opts = list(),
                          fit_opts = list()) {
      # FIXME : to re-enable later
      # # Converting the matrices list to sparse matrix to save space
      # self$A <- lapply(netlist, Matrix::Matrix, sparse = TRUE)
      self$A <- netlist
      
      # Computing the number of rows and cols
      self$nr <- vapply(self$A, nrow, FUN.VALUE = .1)
      self$nc <- vapply(self$A, ncol, FUN.VALUE = .1)

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
      self$model <- model
      self$fit_sbm <- fit_sbm
      self$free_density <-  free_density
      self$free_mixture <- free_mixture
      self$global_opts <- list(Q1_min = 1L,
                               Q1_max = floor(log(sum(self$nr)))+2,
                               Q2_min = 1L,
                               Q2_max = floor(log(sum(self$nc)))+2,
                               sbm_init = TRUE,
                               spectral_init = TRUE,
                               nb_init = 10L,
                               nb_models = 5L,
                               depth = 3L,
                               plot_details = 1L,
                               max_pass = 10L,
                               verbosity = 0L,
                               nb_cores = 1L)
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

      # Initialising the model_list
      self$model_list <- vector("list", self$global_opts$Q1_max * self$global_opts$Q2_max)
      dim(self$model_list) <- c(self$global_opts$Q1_max, self$global_opts$Q2_max)

      if (self$model == "poisson") {
        self$logfactA <- vapply(
          seq_along(self$A),
          function(m) {
            sum(lfactorial(self$A[[m]]) * (1-self$mask[[m]]), na.rm = TRUE)
          },
          FUN.VALUE = .1)
      }
      if(! is.null(self$fit_sbm))  self$global_opts$sbm_init <- FALSE
      self$fit_opts <- list(approx_pois = FALSE,
                            algo_ve = "fp",
                            minibatch = TRUE,
                            verbosity = 0)
      self$fit_opts <- utils::modifyList(self$fit_opts, fit_opts)
    },

    # Fit a list of SBM if fit_sbm == TRUE
    optimize_lbm = function() {
      if (is.null(self$fit_sbm)) {
        self$fit_sbm <-
          bettermc::mclapply(
            X = seq_along(self$A),
            FUN = function(m) {
              sbm::estimateSimpleSBM(
                model = self$model,
                netMat = as.matrix(self$A[[m]]),
                estimOptions = list(verbosity = 0,
                                    plot = FALSE, nbCores = 1L,
                                    exploreMin = self$global_opts$Q_max))
            },
         mc.cores = self$global_opts$nb_cores,
         mc.share.copy = FALSE,
         mc.silent = TRUE
          )
      }
      self$ICL_sbm <- rep(-Inf, self$global_opts$Q_max)
      for (q in seq(self$global_opts$Q_min, self$global_opts$Q_max)) {
        lapply(seq_along(self$A), function(m) self$fit_sbm[[m]]$setModel(index = q))
        self$ICL_sbm[q] <- sum(purrr::MAP_dbl(self$fit_sbm, ~ .$ICL))
      }
      if (self$global_opts$plot_details >= 1) {
        plot(seq(self$global_opts$Q_min, self$global_opts$Q_max),
             self$ICL_sbm[self$global_opts$Q_min : self$global_opts$Q_max],
             col = "green", pch = 5,
             xlab = "Q", ylab = "BICL",
             xlim = c(self$global_opts$Q_min-1, self$global_opts$Q_max+1),
             ylim = c(1.1*max(self$ICL_sbm), .8*max(self$ICL_sbm)))
      }
    },


    optimize_from_sbm = function(index, Q, nb_clusters) {
      # browser()
      # for (q in seq(self$global_opts$Q_min, self$global_opts$Q_max)) {
      lapply(seq_along(self$A), function(m) {
        self$fit_sbm[[m]]$setModel(index = Q)
      })
      # The output clustering from SBM is stored
      # for each of the M models in Z_sbm
      Z_sbm <- lapply(
        seq_along(self$fit_sbm[index]),
        function(m) self$fit_sbm[index][[m]]$memberships
      )
      # The below prob is the intra-cluster connection probability of the alpha
      # ie the diagonal of the alpha matrix
      prob <- lapply(seq_along(self$fit_sbm[index]),
                     function(m) diag(self$fit_sbm[index][[m]]$connectParam$mean))
      nb_init <- ifelse(Q == 1, 1, self$global_opts$nb_init)

      models <- lapply(
        seq(nb_init),
        function(it) {
          if(it == 1) {
            # For the first init, we use the Z order given by the SBM
            mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                           mask = self$mask[index],
                                           model = self$model,
                                           net_id = self$net_id[index],
                                           free_density = self$free_density,
                                           free_mixture = self$free_mixture,
                                           Q = Q,
                                           Z = Z_sbm[index],
                                           logfactA = self$logfactA[index],
                                           init_method = "given",
                                           fit_opts = self$fit_opts)
          } else {
            if(it == 2) {
              Z_init <- lapply(
                seq_along(Z_sbm),
                function(m) {
                  # ord contains Q probabilities and is
                  # deterministically ranked from the lowest to the highest
                  # intra-connection probability
                  ord <- order(prob[[m]])
                  # This returns the cluster membership (Z) in this order
                  # and this clustering is put in Z_init
                  ord[match(Z_sbm[[m]], unique(Z_sbm[[m]]))]
                }
              )
              mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                             mask = self$mask[index],
                                             model = self$model,
                                             net_id = self$net_id,
                                             free_density = self$free_density,
                                             free_mixture = self$free_mixture,
                                             Q = Q,
                                             Z = Z_init,
                                             logfactA = self$logfactA,
                                             # Using the previous re-ordering of Z
                                             # the init_method is "given"
                                             init_method = "given",
                                             fit_opts = self$fit_opts)
            } else {
              Z_init <- lapply(
                seq_along(Z_sbm),
                function(m) {
                  # Here the order is ranked by the highest to the lowest prob
                  # but using a sampling (introducing randomness)
                  ord <- sample(seq_along(prob[[m]]),
                                size = length(prob[[m]]), prob = prob[[m]])
                  ord[match(Z_sbm[[m]], unique(Z_sbm[[m]]))]
                }
              )
              mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                             mask = self$mask[index],
                                             model = self$model,
                                             net_id = self$net_id,
                                             free_density = self$free_density,
                                             free_mixture = self$free_mixture,
                                             Q = Q,
                                             Z = Z_init,
                                             logfactA = self$logfactA,
                                             init_method = "given",
                                             fit_opts = self$fit_opts)
            }
          }
          mypopbm$optimize()
          return(mypopbm)
        }
      )
      best_models <- self$choose_models(models = models, Q = Q)
      self$model_list[[1]][[Q]] <- best_models
      self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
      self$ICL[Q] <- best_models[[1]]$MAP$ICL
      self$BICL[Q] <- best_models[[1]]$BICL
      # }
    },


    optimize_spectral = function(index, Q, nb_clusters) {
      lapply(
        X = seq(self$global_opts$nb_init),
        FUN = function(it) {
          mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                         mask = self$mask[index],
                                         model = self$model,
                                         net_id = self$net_id[index],
                                         free_density = self$free_density,
                                         free_mixture = self$free_mixture,
                                         Q = Q,
                                         logfactA = self$logfactA,
                                         init_method = "spectral",
                                         fit_opts = utils::modifyList(self$fit_opts, list(max_iter = 10L)))
          mypopbm$optimize()
          return(mypopbm)
        }
      )
    },

    optimize_init = function(index, Z, Q, nb_clusters, Cpi = NULL, Calpha = NULL) {
      mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                     mask = self$mask[index],
                                     Z = Z,
                                     model = self$model,
                                     net_id = self$net_id[index],
                                     free_density = self$free_density,
                                     free_mixture = self$free_mixture,
                                     Q = Q,
                                     logfactA = self$logfactA,
                                     init_method = "given",
                                     Cpi = Cpi,
                                     Calpha = Calpha,
                                     fit_opts = self$fit_opts)
      mypopbm$optimize()
      return(mypopbm)
    },

    optimize_from_zinit = function(index, Q, nb_clusters) {
      models <- lapply(self$Z_init[[Q]],
                       function(Z) self$optimize_init(index, Z, Q, nb_clusters))
      best_models <- self$choose_models(models = models, Q = Q)
      self$model_list[[1]][[Q]] <- best_models
      self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
      self$ICL[Q] <- best_models[[1]]$MAP$ICL
      self$BICL[Q] <- best_models[[1]]$BICL
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

      self$separated_inits <- vector("list", 4) # The first coordinate is Q1, the second is Q2
      dim(self$separated_inits) <- c(2,2)

      if (self$global_opts$verbosity >= 3) {
        cat("=== Beginning Burn in ===\n")
      }
        

      # Init for Q = (1,2)
      if (self$global_opts$verbosity >= 4) {
        cat("Fitting ", self$M, " networks for Q = (", toString(c(1, 2)), ")\n")
      }

      self$separated_inits[[1,2]] <- lapply(
        seq.int(self$M),
        function(m) {
          fitBipartiteSBMPop$new(
            A = list(self$A[[m]]), 
            Q = c(1, 2),
            free_mixture = self$free_mixture,
            free_density = self$free_mixture,
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
            A = list(self$A[[m]]), Q = c(2, 1),
            free_mixture = self$free_mixture,
            free_density = self$free_mixture,
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

      # LATER
      # The function fit M LBM for Q = (1,2) and Q = (2,1) from blockmodels
      # TODO : implement the M LBM from blockmodels

      # Here we match the clusters from the M fit objects
      # DONE : implement the matching
      # By using the order of the marginal laws
      if (self$global_opts$verbosity >= 4) {
        cat("\nBeginning to match results for the Separated LBMs.")
      }

      for (m in seq.int(M)) {
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

      # Here we combine the networks to fit a
      # fitBipartite object on the M networks
      M_clusterings_1_2 <- lapply(
        seq.int(M),
        function(m) {
          # We add [[1]] after Z because we fitted
          # only one network with the objects stored
          # where the class is supposed to store more
          self$separated_inits[[1,2]][[m]]$Z[[1]]
        }
      )

      M_clusterings_2_1 <- lapply(
        seq.int(M),
        function(m) {
          # We add [[1]] after Z because we fitted 
          # only one network with the objects stored
          # where the class is supposed to store more
          self$separated_inits[[2,1]][[m]]$Z[[1]]
        }
      )

      self$separated_inits[[1,2]] <- fitBipartiteSBMPop$new(
        A = self$A, Q = c(1, 2),
        free_mixture = self$free_mixture,
        free_density = self$free_mixture,
        Z = M_clusterings_1_2,
        init_method = "given",
        fit_opts = self$fit_opts
      )

      self$separated_inits[[2,1]] <- fitBipartiteSBMPop$new(
        A = self$A, Q = c(2, 1),
        free_mixture = self$free_mixture,
        free_density = self$free_mixture,
        Z = M_clusterings_2_1,
        init_method = "given",
        fit_opts = self$fit_opts
      )

      if (self$global_opts$verbosity >= 4) {
        cat("\nFitting the combined colLBMs.")
      }
      # Here we fit the models
      lapply(seq.int(2),
      function(index){
        # The index used here are a little trick to allow the use
        # of bettermc::mclapply()
        # TODO : parallelize ?
        if (self$global_opts$verbosity >= 4){
          cat(
            "\nFitting the ", self$M, " networks for Q = (",
            toString(c(index, 3 - index)), ")."
          )
        }
        
        self$separated_inits[[index,3 - index]]$optimize()
      })

      if (self$global_opts$verbosity >= 4) {
        cat("\nFinished fitting the colLBM.")
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

      # We parallelize the search from the two points (1,2) / (2,1)
      # and we go looking for the mode with a greedy approach
      # Visiting each of the neighbors
      # greedy_exploration(start)

      

      if(self$global_opts$verbosity >=3) {
        cat(
          "\n==== Finished Burn in",
          " for networks ", self$net_id, " in ", 
          format(Sys.time() - start_time, digits = 3),
          "s ===\n"
        )
        # cat("vbound : ", round(self$vbound), "\n")
        # cat("ICL    : ", round(self$ICL), "\n")
        # cat("BICL   : ", round(self$BICL), "\n")
      }
    },


    forward_pass = function(Q_min = self$global_opts$Q_min,
                            Q_max = self$global_opts$Q_max,
                            index = seq(self$M), nb_clusters = 1L) {
      #browser()
      # if (length(self$BICL) != length(rep(-Inf, Q_max))) {
      #   cat("forward")
      #   cat(self$BICL)
      #   cat(rep(-Inf, Q_max), "\n")
      # }
      old_icl <- pmax(self$BICL, rep(-Inf, self$global_opts$Q_max))
      max_icl <- rep(-Inf, length(old_icl))
      #     max_icl[Q_min] <- self$BICL[Q_min]
      counter <- 0
      #      nb_pass <- 0
      Q <- Q_min +1
      while(Q <= Q_max & counter < self$global_opts$depth) {
        if (is.null(self$model_list[[1]][[Q-1]])) {
          list_popbm <- self$optimize_spectral(index, Q-1, 1L)
          if(! is.null(self$fit_sbm[[Q-1]])) {
            list_popbm <- c(list_popbm, self$optimize_from_sbm(index, Q-1, 1L))
          }
          best_models <- self$choose_models(models = list_popbm, Q = Q-1,
                                            index = index, nb_clusters = nb_clusters)
          best_models <- self$choose_models(models = list_popbm, Q = Q,
                                            index = index, nb_clusters = nb_clusters)
          self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
          self$ICL[Q] <- best_models[[1]]$MAP$ICL
          self$BICL[Q] <- best_models[[1]]$BICL
          max_icl[Q-1] <- best_models[[1]]$BICL

        }
        #        list_popbm <- list()
        model_list <- self$model_list[[1]][[Q-1]]
        #browser()
        list_Zinit <- lapply(
          X = model_list,
          FUN = function(fit) #for (fit in self$model_list[[1]][[Q-1]])
          {
            if (fit$counter_split < 3) {
              fit$counter_split <- fit$counter_split + 1
              Z_init <- purrr::MAP(seq_along(index),
                                   ~ split_clust(as.matrix(fit$A[[.]]), fit$Z[[.]],Q-1))
              Z_init <- purrr::transpose(Z_init)
              return(Z_init)
            } else {
              return (NULL)
            }
          }
        )
        list_Zinit <- Filter(Negate(is.null), list_Zinit)
        list_Zinit <- purrr::flatten(list_Zinit)
        # list_popbm <- parallel::mclapply(
        #   X = model_list,
        #   FUN = function(fit) #for (fit in self$model_list[[1]][[Q-1]])
        #   {
        #     if (fit$counter_split < 3) {
        #       fit$counter_split <- fit$counter_split + 1
        #       Z_init <- purrr::MAP(seq_along(index),
        #                            ~ split_clust(fit$A[[.]], fit$Z[[.]],Q-1))
        #       Z_init <- purrr::transpose(Z_init)
        # k <- min(self$global_opts$nb_cores, ceiling(length(list_Zinit/3)))
    #     browser()
        fold <- split(seq_along(list_Zinit),
                      rep(1:ceiling(length(list_Zinit)/max(Q, 10)),
                          each=max(Q, 10))[1:length(list_Zinit)])
        list_popbm <- bettermc::mclapply(
          seq_along(fold),
          function(x) {
            lapply(
              seq_along(list_Zinit[fold[[x]]]),
              function(it) {
                #browser()
                # res <- self$optimize_init(index = index,
                #                    Z = Z_init[[it]],
                #                    Q = Q,
                #                    nb_clusters = nb_clusters)
                res <- fitSimpleSBMPop$new(A = self$A,
                                           mask = self$mask,
                                           Z = list_Zinit[[it]],
                                           model = self$model,
                                           net_id = self$net_id,
                                           free_density = self$free_density,
                                           free_mixture = self$free_mixture,
                                           Q = Q,
                                           logfactA = self$logfactA,
                                           init_method = "given",
                                           Cpi = NULL,
                                           Calpha = NULL,
                                           fit_opts = self$fit_opts)
                res$optimize()
                tmp_list <- list(res)
                if (self$free_mixture & self$M >1) {
                  C1 <- vapply(seq(res$M),
                               function(m) res$pi[[m]] > 1/res$n[m],
                               FUN.VALUE = rep(TRUE, Q))
                  dim(C1) <- c(Q, res$M)
                  for (zeros in seq(5)) {
                    Cpi <- vapply(seq(res$M),
                                  function(m) res$pi[[m]] > zeros/100, #res$n[m],
                                  FUN.VALUE = rep(TRUE, Q))
                    dim(Cpi) <- c(Q, res$M)
                    if(any(! Cpi) &  all(rowSums(Cpi) > 0) &
                       (zeros == 1 | any(Cpi != C1))) {
                      tmp_res <- res$clone()
                      tmp_res$Cpi <- Cpi
                      tmp_res$Calpha <-
                        Reduce("+", lapply(seq(res$M),
                                           function(m) tcrossprod(Cpi[,m]))) > 0
                      tmp_res$init_method <- "empty"
                      tmp_res$optimize()
                      tmp_list <- c(tmp_list, tmp_res)
                    }
                  }
                }
                # lapply(seq_along(tmp_list),
                #        function(j) {
                #          tmp_list[[j]]$A <- NULL
                #          tmp_list[[j]]$mask <- NULL
                #        }
                # )
                #gc()
                return(tmp_list)
              }
            )
          }, mc.cores = self$global_opts$nb_cores, mc.share.copy = FALSE#, future.globals = list(
          #  model_list = model_list)#,
          #            future.options = list(seed = TRUE)#, mc.cores = 6
        )
        list_popbm <- unlist(list_popbm)
        if (purrr::is_empty(list_popbm) & old_icl[Q] < old_icl[Q-1]) { # a verifier et ou ou
          counter <- counter + 1
        } else {
          best_models <- self$choose_models(models = list_popbm,
                                            Q = Q,
                                            index = index,
                                            nb_clusters = nb_clusters)
          self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
          self$ICL[Q] <- best_models[[1]]$MAP$ICL
          self$BICL[Q] <- best_models[[1]]$BICL
          self$model_list[[nb_clusters]][[Q]] <- best_models
          # lapply(self$model_list[[nb_clusters]][[Q]],
          #        function(fit) {
          #          fit$A <- self$A
          #          fit$mask <- self$mask
          #        }
          # )
          max_icl[Q] <- best_models[[1]]$BICL
          if(self$global_opts$verbosity >= 3) {
            cat(Q, ": ", self$BICL[Q], " -- ", max_icl[Q] > old_icl[Q], "\t")
          }
          if(max_icl[Q] > max_icl[Q-1])  {
            counter <- 0
          } else {
            #   if (max_icl[Q] < old_icl[Q]) {
            counter <- counter + 1
            #}
          }
          Q <- Q+1
        }
        rm(list_popbm)
        #gc()
      }
      #browser()
#      if (any(max_icl > old_icl + .1)) {
      if (max(max_icl) > max(old_icl) + .1) {
        self$improved$forward <- TRUE
      } else {
        self$improved$forward <- FALSE
      }
      Q - 1
    },


    backward_pass = function(Q_min = self$global_opts$Q_min,
                             Q_max = self$global_opts$Q_max,
                             index = seq(self$M), nb_clusters = 1L) {
      # if (length(self$BICL) != length(rep(-Inf, Q_max))) {
      #   cat("backward")
      #   cat(self$BICL)
      #   cat(rep(-Inf, Q_max), "\n")
      # }
      old_icl <- pmax(self$BICL, rep(-Inf, self$global_opts$Q_max))
      max_icl <- rep(-Inf, length(old_icl))

      #     max_icl[Q_max] <- self$BICL[Q_max]
      counter <- 0
      Q <- Q_max - 1
      while(Q >= Q_min & counter < self$global_opts$depth) {
        list_Zinit<- lapply(#furrr::future_MAP(
          X = self$model_list[[1]][[Q+1]],
          FUN = function(fit) {
            if (fit$counter_merge < 3 & all(Reduce("+", fit$MAP$pi) > 0)) {
              fit$counter_merge <- fit$counter_merge + 1
              Z_init <- purrr::MAP(index, ~ merge_clust(fit$Z[[.]], Q+1))
              Z_init <- purrr::transpose(Z_init)
              return(Z_init)
            } else {
              return (NULL)
            }
          }
        )
        list_Zinit <- Filter(Negate(is.null), list_Zinit)
        list_Zinit <- purrr::flatten(list_Zinit)
 #       ind <- sample()
        fold <- split(seq_along(list_Zinit),
                      rep(1:ceiling(length(list_Zinit)/max(Q, 10)), each=max(Q, 10))[1:length(list_Zinit)])
        list_popbm <- bettermc::mclapply(
          seq_along(fold),
          function(x) {
            models <- lapply(
              seq_along(list_Zinit[fold[[x]]]),
              function(it) {
                #browser()
                # res <- self$optimize_init(index = index,
                #                    Z = Z_init[[it]],
                #                    Q = Q,
                #                    nb_clusters = nb_clusters)
                res <- fitSimpleSBMPop$new(A = self$A,
                                           mask = self$mask,
                                           Z = list_Zinit[[it]],
                                           model = self$model,
                                           net_id = self$net_id,
                                           free_density = self$free_density,
                                           free_mixture = self$free_mixture,
                                           Q = Q,
                                           logfactA = self$logfactA,
                                           init_method = "given",
                                           Cpi = NULL,
                                           Calpha = NULL,
                                           fit_opts = self$fit_opts)
                res$optimize()
                tmp_list <- list(res)
                if (self$free_mixture & self$M >1) {
                  C1 <- vapply(seq(res$M),
                               function(m) res$pi[[m]] > 1/res$n[m],
                               FUN.VALUE = rep(TRUE, Q))
                  dim(C1) <- c(Q, res$M)
                  for (zeros in seq(5)) {
                    Cpi <- vapply(seq(res$M),
                                  function(m) res$pi[[m]] > zeros/100, #res$n[m],
                                  FUN.VALUE = rep(TRUE, Q))
                    dim(Cpi) <- c(Q, res$M)
                    if(any(! Cpi) &  all(rowSums(Cpi) > 0) &
                       (zeros == 1 | any(Cpi != C1))) {
                      tmp_res <- res$clone()
                      tmp_res$Cpi <- Cpi
                      tmp_res$Calpha <-
                        Reduce("+", lapply(seq(res$M),
                                           function(m) tcrossprod(Cpi[,m]))) > 0
                      tmp_res$init_method <- "empty"
                      tmp_res$optimize()
                      tmp_list <- c(tmp_list, tmp_res)
                    }
                  }
                }
                # lapply(seq_along(tmp_list),
                #        function(j) {
                #          tmp_list[[j]]$A <- NULL
                #          tmp_list[[j]]$mask <- NULL
                #        }
                # )
                #gc()
                return(tmp_list)
              }
            )
          #  browser()
            models <- unlist(models)
            ord_mod <- order(purrr::MAP_dbl(models, ~.$BICL),
                             decreasing = TRUE)
            best_models <- models[ord_mod[1]]
            for(id in ord_mod) {
              if (length(best_models) >= self$global_opts$nb_models) {
                break
              } else {
                ari <- purrr::MAP_dbl(
                  seq_along(best_models),
                  function(m) {
                    sum(purrr::MAP_dbl(seq_along(self$A),
                                       ~ aricode::ARI(best_models[[m]]$Z[[.]],
                                                      models[[id]]$Z[[.]])))
                  }
                )
                if (all(ari < best_models[[1]]$M)) {
                  best_models <- c(best_models, models[[id]])
                }
              }
            }
            return (best_models)
            }, mc.cores = self$global_opts$nb_cores, mc.share.copy = FALSE#, future.globals = list(
          #  model_list = model_list)#,
          #            future.options = list(seed = TRUE)#, mc.cores = 6
        )
        list_popbm <- unlist(list_popbm)
        if (purrr::is_empty(list_popbm)) {
          counter <- counter + 1
        } else {
          best_models <- self$choose_models(models = list_popbm, Q = Q,
                                            index = index, nb_clusters = nb_clusters)
          self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
          self$ICL[Q] <- best_models[[1]]$MAP$ICL
          self$BICL[Q] <- best_models[[1]]$BICL
          self$model_list[[nb_clusters]][[Q]] <- best_models
          # lapply(self$model_list[[nb_clusters]][[Q]],
          #        function(res) {
          #          res$A <- self$A
          #          res$mask <- self$mask
          #        }
          # )
          max_icl[Q] <- best_models[[1]]$BICL
          if(self$global_opts$verbosity >= 3) {
            cat(Q, ": ", self$BICL[Q], " -- ", max_icl[Q] > old_icl[Q], "\t")
          }
          if(max_icl[Q] > max_icl[Q+1])  {
            counter <- 0
          } else {
          #  if (max_icl[Q] <= old_icl[Q]) {
              counter <- counter + 1
          #  }
          }
        }
        rm(list_popbm)
        #gc()
        Q <- Q - 1
      }
      #browser()
#      if (any(max_icl > old_icl + .1)) {
      if (max(max_icl) > max(old_icl) + .1) {
        self$improved$backward <- TRUE
      } else {
        self$improved$backward <- FALSE
      }
      Q + 1
    },


    optimize = function() {
      #browser()
#      future::plan("future::multisession", workers = self$global_opts$nb_cores)
    #  progressr::handlers(global = TRUE)
     # progressr::handlers("progress")

      # The burn_in step computes models without performing split and merge
      self$burn_in()
      improved <- TRUE
      nb_pass <- 0
      Q <- 1

      # TODO ask @Chabert-Liddell why reduce the number of models ?
      self$global_opts$nb_models <- ceiling(self$global_opts$nb_models/2)
      while (improved & nb_pass < self$global_opts$max_pass) {
        if(self$global_opts$verbosity >=2) {
          cat("==== Starting pass number ", nb_pass + 1,
              " for networks ", self$net_id, " ===\n")
          cat("vbound : ", round(self$vbound), "\n")
          cat("ICL    : ", round(self$ICL), "\n")
          cat("BICL   : ", round(self$BICL), "\n")
        }
        Q <- self$forward_pass(Q_min = max(Q, self$global_opts$Q_min))
        Q <- self$backward_pass(Q_max = min(Q, self$global_opts$Q_max))
        improved <- self$improved$backward | self$improved$forward
        nb_pass <- nb_pass + 1
      }
      self$best_fit <- self$model_list[[1]][[which.max(self$BICL)]][[1]]
      self$model_list[[1]] <- lapply(seq_along(self$model_list[[1]]),
                                     function (Q) self$model_list[[1]][[Q]][1])
    },

    forward_backward = function() {
      #browser()
      counter <- 0
      nb_pass <- 0
      max_icl <- rep(-Inf, self$global_opts$Q_max)
      global_counter <- TRUE
      Q <- self$global_opts$Q_min + 1
      # Forward step
      while (global_counter & nb_pass < self$global_opts$max_pass) {
        global_counter <- FALSE
        while (Q <= self$global_opts$Q_max) {
          list_popbm <- list()
          if (is.null(self$model_list[[1]][[Q-1]])) {
            list_popbm <- self$optimize_spectral(seq(self$M), Q-1, 1L)
            if(! is.null(self$fit_sbm[[Q-1]])) {
              list_popbm <- c(list_popbm, self$optimize_from_sbm(seq(self$M), Q-1, 1L))
            }

          }
          for (fit in self$model_list[[1]][[Q-1]]) {
            Z_init <- purrr::MAP(seq_along(A), ~ split_clust(A[[.]], fit$Z[[.]],Q-1))
            Z_init <- purrr::transpose(Z_init)
            list_res <- lapply(
              seq_along(Z_init),
              function(it) {
                mypopbm <- fitSimpleSBMPop$new(A = A, Q = Q, Z = Z_init[[it]],
                                               logfactA = self$logfactA,
                                               free_density = free_density,
                                               free_mixture = self$free_mixture,
                                               model = model,
                                               fit_opts = self$fit_opts,
                                               init_method = "given")
                mypopbm$optimize()
                return(mypopbm)
              }
            )
            list_popbm <- c(list_popbm, list_res)
          }

          list_spec <-  lapply(
            seq(5),
            function(it) {
              mypopbm <- fitSimpleSBMPop$new(A = A, Q = Q,
                                             logfactA = self$logfactA,
                                             free_density = free_density,
                                             free_mixture = self$free_mixture,
                                             model = model,
                                             fit_opts = self$fit_opts,
                                             init_method = "spectral")
              mypopbm$optimize()
              return(mypopbm)
            }
          )
          list_popbm <- c(list_popbm, list_res, list_spec, best_models[[Q]])
          ord_mod <- order(purrr::MAP_dbl(list_popbm, ~max(.$MAP$ICL, .$BICL)), decreasing = TRUE)
          if(max_icl[Q] < list_popbm[[ord_mod[1]]]$MAP$ICL) global_counter <- TRUE
          max_icl[Q] <- list_popbm[[ord_mod[1]]]$MAP$ICL
          if (max_icl[Q] <= max_icl[Q-1]) counter <- counter + 1
          best_models[[Q]] <- list_popbm[ord_mod[1]]
          # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), ~.$MAP$ICL))
          # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), "BICL"), col = "red")
          for(id in ord_mod) {
            if (length(best_models[[Q]]) < top_models) {
              ari <- purrr::MAP_dbl(
                seq_along(best_models[[Q]]),
                function(m) {
                  sum(purrr::MAP_dbl(seq_along(A),
                                     ~ aricode::ARI(best_models[[Q]][[m]]$Z[[.]],
                                                    list_popbm[[id]]$Z[[.]])))
                }
              )
              if (all(ari < best_models[[Q]][[1]]$M)) {
                best_models[[Q]] <- c(best_models[[Q]], list_popbm[[id]])
                # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), ~.$MAP$ICL))
                # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), "BICL"), col = "red")
              }
            }
          }
          if (counter >= 2) break
          Q <- Q+1
        }
        counter <- 0
        Q <- Q-1

        # Backward step
        while (Q >= max(Q_min,2)) {
          list_popbm <- list()
          for (fit in best_models[[Q+1]]) {
            Z_init <- purrr::MAP(seq_along(A), ~ merge_clust(fit$Z[[.]],Q+1))
            Z_init <- purrr::transpose(Z_init)
            list_res <- lapply(
              seq_along(Z_init),
              function(it) {
                mypopbm <- fitSimpleSBMPop$new(A = A, Q = Q, Z = Z_init[[it]],
                                               logfactA = self$logfactA,
                                               free_density = free_density,
                                               free_mixture = self$free_mixture,
                                               model = model,
                                               fit_opts = self$fit_opts,
                                               init_method = "given",
                                               approx_pois = approx_pois)
                mypopbm$optimize()
                return(mypopbm)
              }
            )
            list_popbm <- c(list_popbm, list_res)
            # points(purrr::MAP_dbl(unlist(list_res), "Q"), purrr::MAP_dbl(unlist(list_res), ~.$MAP$ICL))
            # points(purrr::MAP_dbl(unlist(list_res), "Q"), purrr::MAP_dbl(unlist(list_res), "BICL"), col = "red")
          }
          list_popbm <- c(list_popbm, best_models[[Q]])
          ord_mod <- order(purrr::MAP_dbl(list_popbm, ~max(.$MAP$ICL, .$BICL)), decreasing = TRUE)
          if(max_icl[Q] < list_popbm[[ord_mod[1]]]$MAP$ICL) global_counter <- TRUE

          best_models[[Q]] <- list_popbm[ord_mod[1]]
          max_icl[Q] <- list_popbm[[ord_mod[1]]]$MAP$ICL
          if (max_icl[Q] <= max_icl[Q+1]) counter <- counter + 1
          # points(purrr::MAP_dbl(unlist(best_models), "Q"),
          #        purrr::MAP_dbl(unlist(best_models), ~.$MAP$ICL))
          # points(purrr::MAP_dbl(unlist(best_models), "Q"),
          #        purrr::MAP_dbl(unlist(best_models), "BICL"), col = "red")
          for(id in ord_mod) {
            if (length(best_models[[Q]]) < top_models) {
              ari <- purrr::MAP_dbl(
                seq_along(best_models[[Q]]),
                function(m) {
                  sum(purrr::MAP_dbl(seq_along(A),
                                     ~ aricode::ARI(best_models[[Q]][[m]]$Z[[.]],
                                                    list_popbm[[id]]$Z[[.]])))
                }
              )
              if (all(ari < best_models[[Q]][[1]]$M)) {
                best_models[[Q]] <- c(best_models[[Q]], list_popbm[[id]])
                # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), ~.$MAP$ICL))
                # points(purrr::MAP_dbl(unlist(best_models), "Q"), purrr::MAP_dbl(unlist(best_models), "BICL"), col = "red")
              }
            }
          }
          if (counter >= 2) break
          Q <- Q-1
        }
        Q <- Q+1
      }
      best_id <- which.max(purrr::MAP_dbl(unlist(best_models), ~ .$MAP$ICL))
      unlist(best_models)[[best_id]]
    },

    choose_models = function(models, Q, index = seq(self$M), nb_clusters = 1L) {
      # The provided models are ordered by their BICL in a decreasing order
      ord_mod <- order(purrr::MAP_dbl(models, ~ .$BICL),#~max(.$MAP$ICL, .$BICL)),
                       decreasing = TRUE)

      # If the model_list of the object contains at list Q entries
      # it is appended to the models being processed
      if (length(self$model_list[[nb_clusters]]) >= Q) {
        models <- c(self$model_list[[nb_clusters]][[Q]], models)
      }
      # The models are ordered once again
      # TODO ask @Chabert-Liddell can't it be reduced to just this one step of reordering?
      ord_mod <- order(purrr::MAP_dbl(models, ~.$BICL),
                       decreasing = TRUE)

      # best_models is initialized with the first model of the ord_mod id list
      # ie the one with max BICL
      best_models <- models[ord_mod[1]]
      for(id in ord_mod) {
        # We process the models by their id given by the ord_mod
        if (length(best_models) >= self$global_opts$nb_models) {
          # If we've added the wanted number of models to keep, we exit the loop
          break
        } else {
          # ari is the vector of the model being processed
          # versus all the previously selected best_models
          ari <- purrr::MAP_dbl(
            # This run for each of the best_models
            seq_along(best_models),
            function(m) {
              # Here we sum all the ari for each of the networks clustering
              sum(purrr::MAP_dbl(
                seq_along(self$A),
                ~ aricode::ARI(
                  best_models[[m]]$Z[[.]],
                  models[[id]]$Z[[.]]
                )
              ))
            }
          )
          # If the model has all of his ari less than the number of networks
          # ie the clustering isn't perfect (1 of ARI * M, would be perfect)
          # then the model is  added to the list of best_models
          if (all(ari < best_models[[1]]$M)) {
            best_models <- c(best_models, models[[id]])
          }
        }
      }

      # After having selected the best_models we plot their points
      # x being Q the number of clusters and y being their BICL
      if (self$global_opts$plot_details >= 1) {
        points(purrr::MAP_dbl(unlist(best_models), "Q"),
               purrr::MAP_dbl(unlist(best_models), ~.$BICL))
      }
      return(best_models)
    },


    show = function(type = "Fitted Collection of Bipartite SBM") {
      cat(type, "--", self$model, "variant for", self$M, "networks \n")
      cat("=====================================================================\n")
      cat("net_id = (", self$net_id, ")\n")
      cat(
        "Dimension = (", self$n, ") - (",
        toString(self$best_fit$Q), ") blocks.\n"
      )
      cat("BICL = ", self$best_fit$BICL, " -- #Empty blocks : ", sum(!self$best_fit$Cpi), " \n")
      cat("=====================================================================")
    },

    print = function() self$show(),


    #' Find the points that need to be initialized
    #'
    #' @description
    #' This method look in our models to see if the points for this model exist
    #' @noMd
    #' @noRd
    #' @param center is a vector of the Q1 and Q2 coordinates
    #' in the form c(Q1, Q2)
    #' @param depth is how far away from the center
    #' the function should be applied in a grid style
    #' going from center - (depth,depth) to center + (depth, depth)
    #' @return a list of c(Q1, Q2) for each missing point
  find_missing_points = function(center, depth){
    missing_points <- list()

    model_list <- self$model_list[1]

    center_x  <- center[1]
    center_y <- center[2]

    max_wanted_Q1 <- center_x + depth
    max_wanted_Q2 <- center_y + depth


    # FIXME : take into account that we need the min values
    min_wanted_Q1 <- center_x - depth
    min_wanted_Q2 <- center_y - depth

    missing_points <- append(
      # Here we find the missing points for which at least (Q1,1) is defined
      lapply(
        seq.int(length(model_list)),
        function(q1) {
          # If there are enough block this would be negative
          # so set it to 0 instead, meaning we need
          # 0 more points (Q1,Q2) for this Q1 value
          ifelse(max_wanted_Q2 - length(model_list[[q1]]) >= 0,
            max_wanted_Q2 - length(model_list[[q1]]),
            0
          )
        }
      ),
      # Here we find the missing points for which no point is defined
      ifelse(max_wanted_Q1 > length(model_list),
        # If there we want more Q1 values than there is already
        # we add the wanted_Q2 for each of the missing Q1 values
        rep(max_wanted_Q2, max_wanted_Q1 - length(model_list)),
        # Otherwise we don't need to add anything
        NULL
      )
    )

    return(missing_points)
    },

    #' The moving window application
    #' 
    #' @description
    #' This method is a moving windows 
    #' over the Q1xQ2 space for the number of clusters
    #' @noMd
    #' @noRd
    #' @param center is a vector of the Q1 and Q2 coordinates 
    #' in the form c(Q1, Q2)
    #' @param depth is how far away from the center 
    #' the function should be applied in a grid style 
    #' going from center - (depth,depth) to center + (depth, depth)
    #' @return nothing; but updates the object by adding new models
    moving_window = function(center, depth = 1) {
    # Each split & merge can be parallelized from points
    # But need to be finished when comparing the BICL

    # We loop until there are no missing points in the window

      # Here we compute the missing points
      # TODO : implement the finding of the missing points
      missing_points <- find_missing_points(center, depth)

        # FIXME : Kinda forward_pass
        # Now we compute the possible splits and fit the collection on those points
        # TODO : implement the computation of the splits
        # TODO after : implement the fitting of those splits
        # TODO : select the best one

    # Once all the window is filled 
    # we go to (centerQ1 + depth, centerQ2 + depth) 
    # and go down by merging cluster
    # FIXME : Kinda backward_pass

    }
  )
)
#' An R6 Class object, a collection of model for population of sbm netowrks
#'
#' @export

bmpop <- R6::R6Class(
  "bmpop",
  #
  public = list(
    n = NULL,
    A = NULL,
    M = NULL,
    mask = NULL,
    directed = NULL,
    model = NULL,
    net_id = NULL,
    model_list = NULL,
    global_opts = NULL,
    fit_opts = NULL,
    fit_sbm = NULL,
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


    initialize = function(netlist = NULL,
                          net_id = NULL,
                          directed = NULL,
                          model = "bernoulli",
                          free_density = FALSE,
                          free_mixture = FALSE,
                          fit_sbm = NULL,
                          Z_init = NULL,
                          global_opts = list(),
                          fit_opts = list()) {
      self$A <- lapply(netlist, Matrix::Matrix, sparse=TRUE)
      self$n <- vapply(self$A, nrow, FUN.VALUE = .1)
      self$M <- length(self$A)
      self$mask <- lapply(
        seq_along(self$A),
        function(m) {
          mask <- Matrix::Matrix(diag(1, self$n[m]), sparse = TRUE)
          if (sum(is.na(self$A[[m]])) > 0) {
            mask[is.na(self$A[[m]])] <- 1
          }
          mask
        })
      if (is.null(directed)) {
        self$directed = isSymmetric.matrix(netlist[[1]])
      } else {
        self$directed = directed
      }
      if(is.null(net_id)) {
        self$net_id <- seq(self$M)
      } else {
        self$net_id <- net_id
      }
      self$Z_init <- Z_init
      self$model <- model
      self$fit_sbm <- fit_sbm
      self$free_density <-  free_density
      self$free_mixture <- free_mixture
      self$global_opts <- list(Q_min = 1L,
                               Q_max = floor(log(sum(self$n)))+2,
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
      self$vbound <- rep(-Inf, self$global_opts$Q_max)
      self$ICL <- rep(-Inf, self$global_opts$Q_max)
      self$BICL <- rep(-Inf, self$global_opts$Q_max)
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


    optimize_sbm = function() {
      ## Need to change ICL computation for Z^map
      if (is.null(self$fit_sbm)) {
      #  p <- progressr::progressor(along = self$A)
        self$fit_sbm <-
          bettermc::mclapply(
            X = seq_along(self$A),
            FUN = function(m) {
         #     p(sprintf("m=%g", m))
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
        self$ICL_sbm[q] <- sum(purrr::map_dbl(self$fit_sbm, ~ .$ICL))
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
      #browser()
      # for (q in seq(self$global_opts$Q_min, self$global_opts$Q_max)) {
      lapply(seq_along(self$A), function(m) self$fit_sbm[[m]]$setModel(index = Q))
      Z_sbm <- lapply(seq_along(self$fit_sbm[index]),
                      function(m) self$fit_sbm[index][[m]]$memberships)
      prob <- lapply(seq_along(self$fit_sbm[index]),
                     function(m) diag(self$fit_sbm[index][[m]]$connectParam$mean))
      nb_init <- ifelse(Q == 1, 1, self$global_opts$nb_init)

      models <- lapply(
        seq(nb_init),
        function(it) {
          if(it == 1) {
            mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                           mask = self$mask[index],
                                           model = self$model,
                                           net_id = self$net_id,
                                           directed = self$directed,
                                           free_density = self$free_density,
                                           free_mixture = self$free_mixture,
                                           Q = Q,
                                           Z = Z_sbm,
                                           logfactA = self$logfactA,
                                           init_method = "given",
                                           fit_opts = self$fit_opts)
          } else {
            if(it == 2) {
              Z_init <- lapply(
                seq_along(Z_sbm),
                function(m) {
                  ord <- order(prob[[m]])
                  ord[match(Z_sbm[[m]], unique(Z_sbm[[m]]))]
                }
              )
              mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                             mask = self$mask[index],
                                             model = self$model,
                                             net_id = self$net_id,
                                             directed = self$directed,
                                             free_density = self$free_density,
                                             free_mixture = self$free_mixture,
                                             Q = Q,
                                             Z = Z_init,
                                             logfactA = self$logfactA,
                                             init_method = "given",
                                             fit_opts = self$fit_opts)
            } else {
              Z_init <- lapply(
                seq_along(Z_sbm),
                function(m) {
                  ord <- sample(seq_along(prob[[m]]),
                                size = length(prob[[m]]), prob = prob[[m]])
                  ord[match(Z_sbm[[m]], unique(Z_sbm[[m]]))]
                }
              )
              mypopbm <- fitSimpleSBMPop$new(A = self$A[index],
                                             mask = self$mask[index],
                                             model = self$model,
                                             net_id = self$net_id,
                                             directed = self$directed,
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
      self$ICL[Q] <- best_models[[1]]$map$ICL
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
                                         directed = self$directed,
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
                                     directed = self$directed,
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
      self$ICL[Q] <- best_models[[1]]$map$ICL
      self$BICL[Q] <- best_models[[1]]$BICL
    },

    burn_in = function() {
      # browser()
      if (self$global_opts$sbm_init | ! is.null(self$fit_sbm)) {
        if(self$global_opts$verbosity >=2) {
          cat("Starting optimization of ", self$M, "SBMs")}
        self$optimize_sbm()
      }
      if(! is.null(self$fit_sbm)) {
        if(self$global_opts$verbosity >=2) {
          cat("Starting initialization from SBMs fit. \n")}
        purrr::map(
          .x = seq(self$global_opts$Q_min, self$global_opts$Q_max),
          .f = function(Q) {
            self$optimize_from_sbm(index = seq(self$M),
                                   Q = Q, nb_clusters = 1)}
        )
      }
      if (! is.null(self$Z_init)) {
        if(self$global_opts$verbosity >=2) {
          cat("Starting initialization from given clustering. \n")}
        purrr::map(
          .x = seq(self$global_opts$Q_min, self$global_opts$Q_max),
          .f = function(Q) {
            self$optimize_from_zinit(index = seq(self$M),
                                     Q = Q, nb_clusters = 1)})
      }
      if(self$global_opts$spectral_init) {
        if(self$global_opts$verbosity >=2) {
          cat("Starting initialization from spectral clustering.\n")}
        bettermc::mclapply(
          X = seq(self$global_opts$Q_min, self$global_opts$Q_max),
          FUN = function(Q) {
            models <- self$optimize_spectral(index = seq(self$M),
                                             Q = Q, nb_clusters = 1)
            best_models <- self$choose_models(models = models, Q = Q)
            self$model_list[[1]][[Q]] <- best_models
            self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
            self$ICL[Q] <- best_models[[1]]$map$ICL
            self$BICL[Q] <- best_models[[1]]$BICL
            rm(models)
          }, mc.cores = self$global_opts$nb_cores, mc.share.copy = FALSE
        )
        # voir pour l'init spectral.
      }
      if(self$global_opts$verbosity >=3) {
        cat("==== Finish Burn in",
            " for networks ", self$net_id, " ===\n")
        cat("vbound : ", round(self$vbound), "\n")
        cat("ICL    : ", round(self$ICL), "\n")
        cat("BICL   : ", round(self$BICL), "\n")
      }
      # self$forward_pass()
      # self$backward_pass()
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
          self$ICL[Q] <- best_models[[1]]$map$ICL
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
              Z_init <- purrr::map(seq_along(index),
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
        #       Z_init <- purrr::map(seq_along(index),
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
                                           directed = self$directed,
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
        if (purrr::is_empty(list_popbm) & old_icl[Q] < old_icl[Q-1]) {
          counter <- counter + 1
        } else {
          best_models <- self$choose_models(models = list_popbm,
                                            Q = Q,
                                            index = index,
                                            nb_clusters = nb_clusters)
          self$vbound[Q] <- rev(best_models[[1]]$vbound)[1]
          self$ICL[Q] <- best_models[[1]]$map$ICL
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
        list_Zinit<- lapply(#furrr::future_map(
          X = self$model_list[[1]][[Q+1]],
          FUN = function(fit) {
            if (fit$counter_merge < 3 & all(Reduce("+", fit$map$pi) > 0)) {
              fit$counter_merge <- fit$counter_merge + 1
              Z_init <- purrr::map(index, ~ merge_clust(fit$Z[[.]], Q+1))
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
                                           directed = self$directed,
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
            ord_mod <- order(purrr::map_dbl(models, ~.$BICL),
                             decreasing = TRUE)
            best_models <- models[ord_mod[1]]
            for(id in ord_mod) {
              if (length(best_models) >= self$global_opts$nb_models) {
                break
              } else {
                ari <- purrr::map_dbl(
                  seq_along(best_models),
                  function(m) {
                    sum(purrr::map_dbl(seq_along(self$A),
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
          self$ICL[Q] <- best_models[[1]]$map$ICL
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
      self$burn_in()
      improved <- TRUE
      nb_pass <- 0
      Q <- 1
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
            Z_init <- purrr::map(seq_along(A), ~ split_clust(A[[.]], fit$Z[[.]],Q-1))
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
          ord_mod <- order(purrr::map_dbl(list_popbm, ~max(.$map$ICL, .$BICL)), decreasing = TRUE)
          if(max_icl[Q] < list_popbm[[ord_mod[1]]]$map$ICL) global_counter <- TRUE
          max_icl[Q] <- list_popbm[[ord_mod[1]]]$map$ICL
          if (max_icl[Q] <= max_icl[Q-1]) counter <- counter + 1
          best_models[[Q]] <- list_popbm[ord_mod[1]]
          # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), ~.$map$ICL))
          # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), "BICL"), col = "red")
          for(id in ord_mod) {
            if (length(best_models[[Q]]) < top_models) {
              ari <- purrr::map_dbl(
                seq_along(best_models[[Q]]),
                function(m) {
                  sum(purrr::map_dbl(seq_along(A),
                                     ~ aricode::ARI(best_models[[Q]][[m]]$Z[[.]],
                                                    list_popbm[[id]]$Z[[.]])))
                }
              )
              if (all(ari < best_models[[Q]][[1]]$M)) {
                best_models[[Q]] <- c(best_models[[Q]], list_popbm[[id]])
                # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), ~.$map$ICL))
                # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), "BICL"), col = "red")
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
            Z_init <- purrr::map(seq_along(A), ~ merge_clust(fit$Z[[.]],Q+1))
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
            # points(purrr::map_dbl(unlist(list_res), "Q"), purrr::map_dbl(unlist(list_res), ~.$map$ICL))
            # points(purrr::map_dbl(unlist(list_res), "Q"), purrr::map_dbl(unlist(list_res), "BICL"), col = "red")
          }
          list_popbm <- c(list_popbm, best_models[[Q]])
          ord_mod <- order(purrr::map_dbl(list_popbm, ~max(.$map$ICL, .$BICL)), decreasing = TRUE)
          if(max_icl[Q] < list_popbm[[ord_mod[1]]]$map$ICL) global_counter <- TRUE

          best_models[[Q]] <- list_popbm[ord_mod[1]]
          max_icl[Q] <- list_popbm[[ord_mod[1]]]$map$ICL
          if (max_icl[Q] <= max_icl[Q+1]) counter <- counter + 1
          # points(purrr::map_dbl(unlist(best_models), "Q"),
          #        purrr::map_dbl(unlist(best_models), ~.$map$ICL))
          # points(purrr::map_dbl(unlist(best_models), "Q"),
          #        purrr::map_dbl(unlist(best_models), "BICL"), col = "red")
          for(id in ord_mod) {
            if (length(best_models[[Q]]) < top_models) {
              ari <- purrr::map_dbl(
                seq_along(best_models[[Q]]),
                function(m) {
                  sum(purrr::map_dbl(seq_along(A),
                                     ~ aricode::ARI(best_models[[Q]][[m]]$Z[[.]],
                                                    list_popbm[[id]]$Z[[.]])))
                }
              )
              if (all(ari < best_models[[Q]][[1]]$M)) {
                best_models[[Q]] <- c(best_models[[Q]], list_popbm[[id]])
                # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), ~.$map$ICL))
                # points(purrr::map_dbl(unlist(best_models), "Q"), purrr::map_dbl(unlist(best_models), "BICL"), col = "red")
              }
            }
          }
          if (counter >= 2) break
          Q <- Q-1
        }
        Q <- Q+1
      }
      best_id <- which.max(purrr::map_dbl(unlist(best_models), ~ .$map$ICL))
      unlist(best_models)[[best_id]]
    },

    choose_models = function(models, Q, index = seq(self$M), nb_clusters = 1L) {
      # browser()
      ord_mod <- order(purrr::map_dbl(models, ~ .$BICL),#~max(.$map$ICL, .$BICL)),
                       decreasing = TRUE)
      #      max_icl <- models[[ord_mod[1]]]$map$ICL
      #      icl_improved <- max_icl > self$ICL[[Q]]
      if ( length(self$model_list[[nb_clusters]]) >= Q) {
        models <- c(self$model_list[[nb_clusters]][[Q]], models)
      }
      ord_mod <- order(purrr::map_dbl(models, ~.$BICL),
                       decreasing = TRUE)
      # self$BICL[Q] <- models[ord_mod[1]][[1]]$BICL
      best_models <- models[ord_mod[1]]
      for(id in ord_mod) {
        if (length(best_models) >= self$global_opts$nb_models) {
          break
        } else {
          ari <- purrr::map_dbl(
            seq_along(best_models),
            function(m) {
              sum(purrr::map_dbl(seq_along(self$A),
                                 ~ aricode::ARI(best_models[[m]]$Z[[.]],
                                                models[[id]]$Z[[.]])))
            }
          )
          if (all(ari < best_models[[1]]$M)) {
            best_models <- c(best_models, models[[id]])
          }
        }
      }
      if (self$global_opts$plot_details >= 1) {
        points(purrr::map_dbl(unlist(best_models), "Q"),
               purrr::map_dbl(unlist(best_models), ~.$BICL))
      }
      return(best_models)
    },


    show = function(type = "Fitted Collection of Simple SBM") {
      cat(type, "--", self$model, "variant for", self$M, "networks \n")
      cat("=====================================================================\n")
      cat("net_id = (", self$net_id, ")\n")
      cat("Dimension = (", self$n, ") - (",
          self$best_fit$Q,  ") blocks.\n")
      cat("BICL = ", self$best_fit$BICL, " -- #Empty blocks : ", sum(!self$best_fit$Cpi), " \n")
      cat("=====================================================================")
    },

    print = function() self$show()
  )
)




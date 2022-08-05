#' Estimate a colSBM on a collection of networks
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "delta",
#' "deltapi".
#' @param net_id A vector of string, the name of the networks.
#' @param directed A boolean, are the networks directed or not.
#' @param model A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param fit_sbm A list of model using the \code{sbm} package. Use to speed up
#' the initialization.
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and the outpur
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#'
#' @return A bmpop object listing a collection of models for the collection.
#' of networks
#' @export
#'
#' @examples
estimate_colSBM <-
  function(    netlist,
               colsbm_model,
               net_id = NULL,
               directed = NULL,
               model = "bernoulli",
               fit_sbm = NULL,
               nb_run = 3L,
               global_opts = list(),
               fit_opts = list(),
               fit_init = NULL) {
    switch (colsbm_model,
            "iid" = {
              free_density <-  FALSE
              free_mixture <-  FALSE
            },
            "pi" = {
              free_density <-  FALSE
              free_mixture <-  TRUE
            },
            "delta" = {
              free_density <-  TRUE
              free_mixture <-  FALSE
            },
            "deltapi" = {
              free_density <-  TRUE
              free_mixture <-  TRUE
            },
            stop("colsbm_model unknown. Must be one of iid, pi, delta or deltapi"))
    if (is.null(global_opts$nb_cores)) {
      global_opts$nb_cores <- 1L
    }
    nb_cores <- global_opts$nb_cores
    if (is.null(global_opts$Q_max)) {
      Q_max <- floor(log(sum(sapply(netlist, function(A) nrow(A))))+2)
    } else {
      Q_max <- global_opts$Q_max
    }
    if (! is.null(fit_init)) {
      my_bmpop <- fit_init
    } else {
      if(is.null(net_id)) {
        net_id <- seq_along(netlist)
      }
      if(is.null(fit_sbm)) {
        fit_sbm <-
          bettermc::mclapply(
            X = seq_along(netlist),
            FUN = function(m) {
              #     p(sprintf("m=%g", m))
              sbm::estimateSimpleSBM(
                model = model,
                netMat = netlist[[m]],
                estimOptions = list(verbosity = 0,
                                    plot = FALSE, nbCores = 1L,
                                    exploreMin = Q_max))
            }, mc.cores = nb_cores, mc.silent = TRUE
          )
      }
      tmp_fits <-
        bettermc::mclapply(
          seq(nb_run),
          function(x) {
            global_opts$nb_cores <- max(1L, floor(global_opts$nb_cores/nb_run))
            tmp_fit <- bmpop$new(netlist = netlist,
                                 net_id = net_id,
                                 directed = directed,
                                 model = model,
                                 free_density = free_density,
                                 free_mixture = free_mixture,
                                 fit_sbm = fit_sbm,
                                 global_opts = global_opts,
                                 fit_opts = fit_opts)
            tmp_fit$optimize()
            return(tmp_fit)
          }, mc.progress = TRUE, mc.cores = min(nb_run,nb_cores)
        )
      my_bmpop <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$BICL,
                                             FUN.VALUE = .1))]]
      my_bmpop$model_list[[1]] <-
        lapply(X = seq_along(my_bmpop$model_list[[1]]),
               FUN = function(q) {
                 tmp_fits[[which.max(vapply(
                   tmp_fits,
                   function(fit) fit$model_list[[1]][[q]][[1]]$BICL,
                   FUN.VALUE = .1))]]$model_list[[1]][[q]]
               }
        )
      my_bmpop$vbound <- vapply(my_bmpop$model_list[[1]],
                             function (fit) rev(fit[[1]]$vbound)[1], FUN.VALUE = .1)
      my_bmpop$ICL <- vapply(my_bmpop$model_list[[1]],
                             function (fit) fit[[1]]$ICL, FUN.VALUE = .1)
      my_bmpop$BICL <- vapply(my_bmpop$model_list[[1]],
                                        function (fit) fit[[1]]$BICL, FUN.VALUE = .1)
      if(my_bmpop$global_opts$verbosity >=1) {
        cat("==== Optimization finished for networks ", my_bmpop$net_id, " ====\n")
        cat("vbound : ", round(my_bmpop$vbound), "\n")
        cat("ICL    : ", round(my_bmpop$ICL), "\n")
        cat("BICL   : ", round(my_bmpop$BICL), "\n")
        cat("Best model for Q = ", which.max(my_bmpop$BICL), "\n")
        if (! is.null(my_bmpop$ICL_sbm)) {
          icl_sbm <- vapply(seq(my_bmpop$M),
                            function(m) max(my_bmpop$fit_sbm[[m]]$storedModels$ICL),
                            FUN.VALUE = .1)
          if (max(my_bmpop$BICL) > sum(icl_sbm)) {
            cat("Joint modelisation preferred over separated one. BICL: " )
          } else {
            cat("Joint modelisation not preferred over separated one. BICL: " )
          }
          cat(max(my_bmpop$BICL), " vs ", sum(icl_sbm))
        }
      }
      rm(tmp_fits)
      gc()
    }
    return(my_bmpop)
  }

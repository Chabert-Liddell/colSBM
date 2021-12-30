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
            }, mc.cores = nb_cores
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
      my_bmpop <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$ICL_clustering,
                                             FUN.VALUE = .1))]]
      my_bmpop$model_list[[1]] <-
        lapply(X = seq_along(my_bmpop$model_list[[1]]),
               FUN = function(q) {
                 tmp_fits[[which.max(vapply(
                   tmp_fits,
                   function(fit) fit$model_list[[1]][[q]][[1]]$ICL_clustering,
                   FUN.VALUE = .1))]]$model_list[[1]][[q]]
               }
        )
      my_bmpop$ICL <- vapply(my_bmpop$model_list[[1]],
                             function (fit) fit[[1]]$ICL, FUN.VALUE = .1)
      my_bmpop$ICL_clustering <- vapply(my_bmpop$model_list[[1]],
                                        function (fit) fit[[1]]$ICL_clustering, FUN.VALUE = .1)
      rm(tmp_fits)
      gc()
    }
    return(my_bmpop)
  }

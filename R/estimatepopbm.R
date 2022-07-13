

#' Title
#'
#' @param netlist
#' @param net_id
#' @param directed
#' @param model
#' @param free_density
#' @param free_mixture
#' @param fit_sbm
#' @param global_opts
#' @param fit_opts
#'
#' @return
#' @export
#'
#' @examples
EstimatePopBM <- function(netlist = NULL,
                          net_id = NULL,
                          directed = NULL,
                          model = "bernoulli",
                          free_density = FALSE,
                          free_mixture = TRUE,
                          fit_sbm = NULL,
                          nb_run = 3L,
                          global_opts = list(),
                          fit_opts = list(),
                          full_inference = FALSE,
                          fit_init = NULL) {

  ## Ajouter une sorte de smoothing par les clusters de reseaux
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



  # list_model_binary <- list(my_bmpop)
  # list_icl_binary <- list(my_bmpop$best_fit$ICL)
  recursive_clustering <- function(fit) {
    if(fit$best_fit$M == 1) return(fit$best_fit)
    dist_bm <- diag(0, fit$best_fit$M)
    for (i in seq(nrow(dist_bm))) {
      for (j in seq(ncol(dist_bm))) {
        if(free_mixture) {
          fit_pi <- fit$best_fit$pi[c(i,j)]
        } else {
          fit_pi <- fit$best_fit$pim[c(i,j)]
        }
        dist_bm[i,j] <-
          dist_bmpop_max(
            pi = fit_pi,
            alpha = fit$best_fit$alpham[c(i,j)],
            delta = fit$best_fit$delta[c(i,j)], weight = "max")
      }
    }
  #    cl <- stats::cutree(stats::hclust(as.dist(dist_bm), method = "ward.D2"), 2)
      if (fit$M >= 3) {
        cl <- cluster::pam(x = sqrt(dist_bm), k = 2, diss = TRUE)$clustering
      } else {
        cl <- stats::cutree(stats::hclust(as.dist(dist_bm), method = "ward.D2"), 2)
      }

    fits <-
      lapply(c(1,2),
             function(k) {
               Z_init <- lapply(seq_along(fit$model_list[[1]]),
                                function(q) {
                                  lapply(seq_along(fit$model_list[[1]][[q]]),
                                         function(j) fit$model_list[[1]][[q]][[j]]$Z[cl == k])
                                })
               if (is.null(global_opts$nb_cores)) {
                 global_opts$nb_cores <- 1L
               } else {
                 nb_cores <- global_opts$nb_cores
               }
               tmp_fits <-
                 bettermc::mclapply(
                   seq(min(sum(cl == k), nb_run)),
                   function(x) {
                     global_opts$nb_cores <- max(1L, floor(global_opts$nb_cores/nb_run))
                     tmp_fit <- bmpop$new(fit$A[cl == k], fit$net_id[cl == k],
                                      directed = directed,
                                      model = model,
                                      free_density = free_density,
                                      free_mixture = free_mixture,
                                      fit_sbm = fit$fit_sbm[cl==k],
                                      Z_init = Z_init,
                                      global_opts = global_opts,
                                      fit_opts = fit_opts)
                     # global_opts = list(nb_models = 10L,
                     #                    nb_init = 100, verbosity = 1))
                     tmp_fit$optimize()
                     return(tmp_fit)
                   }, mc.progress = TRUE, mc.cores = min(nb_run,nb_cores)
                   )
               res <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$ICL_clustering,
                                         FUN.VALUE = .1))]]
               res$model_list[[1]] <-
                 lapply(X = seq_along(res$model_list[[1]]),
                        FUN = function(q) {
                          tmp_fits[[which.max(vapply(
                            tmp_fits,
                            function(fit) fit$model_list[[1]][[q]][[1]]$ICL_clustering,
                            FUN.VALUE = .1))]]$model_list[[1]][[q]]
                        }
                 )
               res$ICL <- vapply(res$model_list[[1]],
                                      function (fit) fit[[1]]$ICL, FUN.VALUE = .1)
               res$ICL_clustering <- vapply(res$model_list[[1]],
                                      function (fit) fit[[1]]$ICL_clustering, FUN.VALUE = .1)
               rm(tmp_fits)
               gc()
               return(res)
             })
    if (full_inference) {
      return(list(fit$best_fit,
                  recursive_clustering(fits[[1]]),
                  recursive_clustering(fits[[2]])))
    } else {
      if (fits[[1]]$best_fit$ICL_clustering + fits[[2]]$best_fit$ICL_clustering >
          fit$best_fit$ICL_clustering) {
        return(list(fit$best_fit,
                    recursive_clustering(fits[[1]]),
                    recursive_clustering(fits[[2]])))
      } else {
        return(fit$best_fit)
      }
    }
    #   return (list(list(model = model, id = model$net_id, ICL = model$best_fit$ICL, Q = model$best_fit$Q),
    #                list(list(model = models[[1]], id = models[[1]]$net_id, ICL = models[[1]]$best_fit$ICL, Q = models[[1]]$best_fit$Q),
    #                     list(model = models[[2]], id = models[[2]]$net_id, ICL = models[[2]]$best_fit$ICL, Q = models[[2]]$best_fit$Q))))
  }

  list_model_binary <- recursive_clustering(my_bmpop)


  #======================= Clustering par empty cluster

  # rk <- rev(rank(my_bmpop$ICL_clustering))
  # net_cl <- purrr::map(rk, ~my_bmpop$model_list[[1]][[.]][[1]]$net_clustering)
  # net_cl <- unique(net_cl)
  # net_cl <- net_cl[1:(min(length(net_cl), 5))]
  # part_id <- list()
  # for(nc in seq_along(net_cl)) {
  #   tmp_id <- lapply(unique(net_cl[[nc]]), function(g) net_id[net_cl[[nc]] == g])
  #   part_id <- c(part_id, tmp_id)
  # }
  # #browser()
  # id_list <- vector("list", length(net_cl))
  # partial_id <- unique(c(list(net_id), part_id))
  # for(nc in seq_along(net_cl)) {
  #   id_list[[nc]] <-
  #     c(id_list[[nc]],
  #       sapply(unique(net_cl[[nc]]),
  #              function(g) {
  #                x <- net_id[net_cl[[nc]] == g]
  #                which(
  #                  unlist(
  #                    lapply(
  #                      seq_along(partial_id),
  #                      function(f) isTRUE(all.equal(x, partial_id[[f]]))
  #                    )))
  #              }
  #       ))
  # }
  # partial_fit <- list(my_bmpop)
  # it <- 2
  # #browser()
  # while (it <= length(partial_id)) {
  #   Z_init <- lapply(
  #     seq_along(my_bmpop$model_list[[1]]),
  #     function(q) {
  #       lapply(seq_along(my_bmpop$model_list[[1]][[q]]),
  #              function(j) my_bmpop$model_list[[1]][[q]][[j]]$Z[partial_id[[it]]])
  #                    }
  #       )
  #   partial_fit[[it]] <-
  #     bmpop$new(netlist = netlist[partial_id[[it]]],
  #               net_id = net_id[partial_id[[it]]],
  #               directed = directed,
  #               model = model,
  #               free_density = free_density,
  #               Z_init = Z_init,
  #               fit_sbm = fit_sbm[partial_id[[it]]],
  #               global_opts = global_opts,
  #               fit_opts = fit_opts)
  #   partial_fit[[it]]$optimize()
  #   tmp_rk <- rev(rank(partial_fit[[it]]$ICL_clustering))
  #   tmp_cl <- purrr::map(tmp_rk,
  #                        ~ partial_fit[[it]]$model_list[[1]][[.]][[1]]$net_clustering)
  #
  #   tmp_cl <- unique(tmp_cl)
  #   tmp_cl <- tmp_cl[1:(min(length(tmp_cl), 5))]
  #   tmp_id <- list()
  #   for (nc in seq_along(tmp_cl)) {
  #     tmp_id <- c(tmp_id,
  #                     lapply(unique(tmp_cl[[nc]]), function(g) partial_id[[it]][tmp_cl[[nc]] == g]))
  #   }
  #   old_id <- length(partial_id)
  #   partial_id <- unique(c(partial_id, tmp_id))
  #   new_id <- length(partial_id)
  #
  #   if (old_id < new_id) {
  #     index <- which(sapply(seq_along(id_list), function(q) any(id_list[[q]] == it)))
  #     id_list <- c(id_list,
  #                  lapply(
  #                    seq_along(id_list[index]),
  #                    function(q) {
  #                      idl <- id_list[index][[q]]
  #                      wh <- numeric(new_id - old_id +1)
  #                      for (id in seq(old_id+1, new_id)) {
  #                        vv <- which(vapply(seq_along(partial_id),
  #                                           function(k)
  #                                             setequal(setdiff(partial_id[[it]],
  #                                                              partial_id[[id]]),
  #                                                      partial_id[[k]]), FUN.VALUE = TRUE))
  #                        wh[id-old_id] <- ifelse(length(vv) > 0, vv, 0)
  #                      }
  #                      sort(unique(c(idl[idl != it], na.omit(wh[wh!=0]), (old_id+1):new_id)))
  #     }))
  #   }
  #   it <- it + 1
  # }
  # return(list(id_list = id_list,
  #             id = partial_id,
  #             fit = partial_fit))

}




best_list <- function(l) {
  if(length(l) == 1)  return(l[[1]])
  if(length(l) == 2) return(l[[2]])
  return(list(best_list(l[[2]]),
              best_list(l[[3]])))
}

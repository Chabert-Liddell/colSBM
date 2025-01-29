#' Partition of a collection of networks based on their common
#' mesoscale structures
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "delta",
#' "deltapi".
#' @param net_id A vector of string, the name of the networks.
#' @param directed A boolean, are the networks directed or not.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param fit_sbm A list of fitted models using the \code{sbm} package.
#' Use to speed up the initialization.
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and the output
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#' @param full_inference The default "FALSE", the algorithm stop once splitting
#' groups of networks does not improve the BICL criterion. If "TRUE", then
#' continue to split groups until a trivial classification of one network per
#' group.
#'
#' @return A list of models for the recursive partition of
#' the collection of networks.
#'
#' @details The best partition could be extract with the function
#' `extract_best_partition()`. The object of the list are FitSimpleSBMPop object,
#' so it is a model for a given number of blocks Q.
#' @export
#'
#' @seealso [colSBM::extract_best_partition()], [colSBM::estimate_colSBM()],
#' \code{\link[colSBM]{fitSimpleSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#'
#' #' # Trivial example with Gnp networks:
#' Net <- lapply(
#'   list(.7, .7, .2, .2),
#'   function(p) {
#'     A <- matrix(0, 15, 15)
#'     A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
#'     A <- A + t(A)
#'   }
#' )
#' \dontrun{
#' cl <- clusterize_networks(Net,
#'   colsbm_model = "iid",
#'   directed = FALSE,
#'   distribution = "bernoulli",
#'   nb_run = 1
#' )
#' }
clusterize_networks <- function(netlist,
                                colsbm_model,
                                net_id = NULL,
                                directed = NULL,
                                distribution = "bernoulli",
                                fit_sbm = NULL,
                                nb_run = 3L,
                                global_opts = list(),
                                fit_opts = list(),
                                fit_init = NULL,
                                full_inference = FALSE) {
  ## Ajouter une sorte de smoothing par les clusters de reseaux
  ##
  ##
  switch(colsbm_model,
    "iid" = {
      free_density <- FALSE
      free_mixture <- FALSE
    },
    "pi" = {
      free_density <- FALSE
      free_mixture <- TRUE
    },
    "delta" = {
      free_density <- TRUE
      free_mixture <- FALSE
    },
    "deltapi" = {
      free_density <- TRUE
      free_mixture <- TRUE
    },
    stop("colsbm_model unknown. Must be one of iid, pi, delta or deltapi")
  )



  if (is.null(global_opts$nb_cores)) {
    global_opts$nb_cores <- 1L
  }
  nb_cores <- global_opts$nb_cores
  if (is.null(global_opts$Q_max)) {
    Q_max <- floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)
  } else {
    Q_max <- global_opts$Q_max
  }

  if (is.null(global_opts$backend)) {
    global_opts$backend <- "parallel"
  }

  if (!is.null(fit_init)) {
    my_bmpop <- fit_init
  } else {
    if (is.null(net_id)) {
      net_id <- seq_along(netlist)
    }
    if (is.null(fit_sbm)) {
      fit_sbm <-
        colsbm_lapply(
          X = seq_along(netlist),
          FUN = function(m) {
            #     p(sprintf("m=%g", m))
            sbm::estimateSimpleSBM(
              model = distribution,
              netMat = netlist[[m]],
              estimOptions = list(
                verbosity = 0,
                plot = FALSE, nbCores = 1L,
                exploreMin = Q_max
              )
            )
          },
          backend = global_opts$backend,
          nb_cores = nb_cores
        )
    }
    tmp_fits <-
      colsbm_lapply(
        seq(nb_run),
        function(x) {
          global_opts$nb_cores <- max(1L, floor(global_opts$nb_cores / nb_run))
          tmp_fit <- bmpop$new(
            netlist = netlist,
            net_id = net_id,
            directed = directed,
            distribution = distribution,
            free_density = free_density,
            free_mixture = free_mixture,
            fit_sbm = fit_sbm,
            global_opts = global_opts,
            fit_opts = fit_opts
          )
          tmp_fit$optimize()
          return(tmp_fit)
        },
        backend = global_opts$backend,
        nb_cores = min(nb_run, nb_cores),
        mc.progress = TRUE # , nb_cores = min(nb_run, nb_cores)
      )
    my_bmpop <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$BICL,
      FUN.VALUE = .1
    ))]]
    my_bmpop$model_list[[1]] <-
      lapply(
        X = seq_along(my_bmpop$model_list[[1]]),
        FUN = function(q) {
          tmp_fits[[which.max(vapply(
            tmp_fits,
            function(fit) fit$model_list[[1]][[q]][[1]]$BICL,
            FUN.VALUE = .1
          ))]]$model_list[[1]][[q]]
        }
      )
    my_bmpop$ICL <- vapply(my_bmpop$model_list[[1]],
      function(fit) fit[[1]]$ICL,
      FUN.VALUE = .1
    )
    my_bmpop$BICL <- vapply(my_bmpop$model_list[[1]],
      function(fit) fit[[1]]$BICL,
      FUN.VALUE = .1
    )
    rm(tmp_fits)
    gc()
  }



  # list_model_binary <- list(my_bmpop)
  # list_icl_binary <- list(my_bmpop$best_fit$ICL)
  recursive_clustering <- function(fit) {
    if (fit$best_fit$M == 1) {
      return(fit$best_fit)
    }
    dist_bm <- diag(0, fit$best_fit$M)
    for (i in seq(nrow(dist_bm))) {
      for (j in seq(ncol(dist_bm))) {
        if (free_mixture) {
          fit_pi <- fit$best_fit$pi[c(i, j)]
        } else {
          fit_pi <- fit$best_fit$pim[c(i, j)]
        }
        dist_bm[i, j] <-
          dist_bmpop_max(
            pi = fit_pi,
            alpha = fit$best_fit$alpham[c(i, j)],
            delta = fit$best_fit$delta[c(i, j)], weight = "max"
          )
      }
    }
    #    cl <- stats::cutree(stats::hclust(as.dist(dist_bm), method = "ward.D2"), 2)
    if (fit$M >= 3) {
      cl <- cluster::pam(x = sqrt(dist_bm), k = 2, diss = TRUE)$clustering
    } else {
      cl <- stats::cutree(stats::hclust(stats::as.dist(dist_bm), method = "ward.D2"), 2)
    }

    fits <-
      lapply(
        c(1, 2),
        function(k) {
          Z_init <- lapply(
            seq_along(fit$model_list[[1]]),
            function(q) {
              lapply(
                seq_along(fit$model_list[[1]][[q]]),
                function(j) fit$model_list[[1]][[q]][[j]]$Z[cl == k]
              )
            }
          )
          if (is.null(global_opts$nb_cores)) {
            global_opts$nb_cores <- 1L
          } else {
            nb_cores <- global_opts$nb_cores
          }
          tmp_fits <-
            colsbm_lapply(
              seq(min(sum(cl == k), nb_run)),
              function(x) {
                global_opts$nb_cores <- max(1L, floor(global_opts$nb_cores / nb_run))
                tmp_fit <- bmpop$new(fit$A[cl == k], fit$net_id[cl == k],
                  directed = directed,
                  distribution = distribution,
                  free_density = free_density,
                  free_mixture = free_mixture,
                  fit_sbm = fit$fit_sbm[cl == k],
                  Z_init = Z_init,
                  global_opts = global_opts,
                  fit_opts = fit_opts
                )
                # global_opts = list(nb_models = 10L,
                #                    nb_init = 100, verbosity = 1))
                tmp_fit$optimize()
                return(tmp_fit)
              },
              nb_cores = min(nb_run, nb_cores),
              backend = global_opts$backend
            )
          res <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$BICL,
            FUN.VALUE = .1
          ))]]
          res$model_list[[1]] <-
            lapply(
              X = seq_along(res$model_list[[1]]),
              FUN = function(q) {
                tmp_fits[[which.max(vapply(
                  tmp_fits,
                  function(fit) fit$model_list[[1]][[q]][[1]]$BICL,
                  FUN.VALUE = .1
                ))]]$model_list[[1]][[q]]
              }
            )
          res$ICL <- vapply(res$model_list[[1]],
            function(fit) fit[[1]]$ICL,
            FUN.VALUE = .1
          )
          res$BICL <- vapply(res$model_list[[1]],
            function(fit) fit[[1]]$BICL,
            FUN.VALUE = .1
          )
          rm(tmp_fits)
          gc()
          return(res)
        }
      )
    if (full_inference) {
      return(list(
        fit$best_fit,
        recursive_clustering(fits[[1]]),
        recursive_clustering(fits[[2]])
      ))
    } else {
      if (fits[[1]]$best_fit$BICL + fits[[2]]$best_fit$BICL >
        fit$best_fit$BICL) {
        return(list(
          fit$best_fit,
          recursive_clustering(fits[[1]]),
          recursive_clustering(fits[[2]])
        ))
      } else {
        return(fit$best_fit)
      }
    }
  }

  list_model_binary <- recursive_clustering(my_bmpop)

  invisible(list_model_binary)
}





#' Extract the best partition from the list of model given by the functions
#' `clusterize_networks()` or `clusterize_bipartite_networks()`.
#'
#' @param l A list of models obtained from the function  `clusterize_networks()`
#' @param unnest A boolean specifying if the returned object should be un-nested
#' (and thus loose exploration clustering structure) or not. Default to TRUE.
#'
#'
#' @return A list of models giving the best partition.
#' @export
#'
#' @examples
#'
#' #' # Trivial example with Gnp networks:
#' Net <- lapply(
#'   list(.7, .7, .2, .2),
#'   function(p) {
#'     A <- matrix(0, 15, 15)
#'     A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
#'     A <- A + t(A)
#'   }
#' )
#' \dontrun{
#' cl <- clusterize_networks(Net,
#'   colsbm_model = "iid",
#'   directed = FALSE,
#'   distribution = "bernoulli",
#'   nb_run = 1
#' )
#' best_partition <- extract_best_partition(cl)
#' }
extract_best_partition <- function(l, unnest = TRUE) {
  if (inherits(l, "fitSimpleSBMPop") || inherits(l, "fitBipartiteSBMPop")) {
    return(l)
  }
  stopifnot("Provided clustering contains incorrect fit objects ! Should be fitSimpleSBMPop or fitBipartiteSBMPop" = (inherits(l[[1]], "fitSimpleSBMPop") || inherits(l[[1]], "fitBipartiteSBMPop")))
  if (length(l) == 1) {
    return(l[[1]])
  }
  if (length(l) == 2) {
    return(l[[2]])
  }
  out <- list(
    extract_best_partition(l[[2]]),
    extract_best_partition(l[[3]])
  )
  if (unnest && is.list(out)) {
    out <- unlist(out)
  }
  return(out)
}

#' Partition of a collection of bipartite networks based on their common
#' mesoscale structures
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colBiSBM to use, one of "iid", "pi", "rho", "pirho",
#' @param net_id A vector of string, the name of the networks.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and the output
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#' @param full_inference The default "FALSE", the algorithm stop once splitting
#' groups of networks does not improve the BICL criterion. If "TRUE", then
#' continue to split groups until a trivial classification of one network per
#' group.
#' @param verbose A boolean, should the function be verbose or not. Default to
#' TRUE.
#'
#' @importFrom future.apply future_lapply
#' @import cli
#' @importFrom utils modifyList
#'
#' @return A list of models for the recursive partition of
#' the collection of networks.
#'
#' @details The best partition could be extract with the function
#' `extract_best_partition()`. The object of the list are fitBipartiteSBMPop
#' object, so it is a model for a given number of blocks Q1, Q2.
#'
#' This functions make call to `estimate_colBiSBM`.
#' @export
#'
#' @seealso [colSBM::extract_best_partition()], [colSBM::estimate_colBiSBM()],
#' \code{\link[colSBM]{fitBipartiteSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#' alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
#' alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
#' first_collection <- generate_bipartite_collection(
#'   nr = 50, nc = 25,
#'   pi = c(0.5, 0.5), rho = c(0.5, 0.5),
#'   alpha = alpha1, M = 2
#' )
#' second_collection <- generate_bipartite_collection(
#'   nr = 50, nc = 25,
#'   pi = c(0.5, 0.5), rho = c(0.5, 0.5),
#'   alpha = alpha2, M = 2
#' )
#'
#' netlist <- append(first_collection, second_collection)
#'
#' \dontrun{
#' cl_separated <- clusterize_bipartite_networks(
#'   netlist = netlist,
#'   colsbm_model = "iid",
#'   global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
#' )
#' }
clusterize_bipartite_networks <- function(netlist,
                                          colsbm_model,
                                          net_id = NULL,
                                          distribution = "bernoulli",
                                          nb_run = 3L,
                                          global_opts = list(),
                                          fit_opts = list(),
                                          fit_init = NULL,
                                          full_inference = FALSE,
                                          verbose = TRUE) {
  # Adding default global_opts
  switch(colsbm_model,
    "iid" = {
      free_mixture_row <- FALSE
      free_mixture_col <- FALSE
    },
    "pi" = {
      free_mixture_row <- TRUE
      free_mixture_col <- FALSE
    },
    "rho" = {
      free_mixture_row <- FALSE
      free_mixture_col <- TRUE
    },
    "pirho" = {
      free_mixture_row <- TRUE
      free_mixture_col <- TRUE
    },
    stop(
      "colsbm_model unknown.",
      " Must be one of iid, pi, rho, pirho, delta or deltapi"
    )
  )
  # Check if a netlist is provided, try to cast it if not
  if (!is.list(netlist)) {
    netlist <- list(netlist)
  }
  # go is used to temporarily store the default global_opts
  go <- list(
    Q1_min = 1L,
    Q2_min = 1L,
    Q1_max = floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2),
    Q2_max = floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2),
    nb_init = 10L,
    nb_models = 5L,
    backend = "future",
    depth = 1L,
    plot_details = 0L,
    max_pass = 10L,
    verbosity = 1L,
    nb_cores = 1L
  )
  go <- utils::modifyList(go, global_opts)
  global_opts <- go
  if (is.null(global_opts$nb_cores)) {
    global_opts$nb_cores <- 1L
  }
  nb_cores <- global_opts$nb_cores
  if (is.null(global_opts$backend)) {
    global_opts$backend <- "parallel"
  }
  if (is.null(global_opts$Q1_max)) {
    Q1_max <- floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)
  } else {
    Q1_max <- global_opts$Q1_max
  }
  if (is.null(global_opts$Q2_max)) {
    Q2_max <- floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2)
  } else {
    Q2_max <- global_opts$Q2_max
  }

  # Fit the initial model on the full collection
  if (verbose) {
    cli::cli_h1("Fitting the full collection")
  }

  my_bisbmpop <- estimate_colBiSBM(
    netlist = netlist,
    colsbm_model = colsbm_model,
    net_id = net_id,
    distribution = distribution,
    nb_run = nb_run,
    global_opts = global_opts,
    fit_opts = fit_opts
  )

  clustering_queue <- list(my_bisbmpop)
  list_model_binary <- list()

  if (verbose) {
    cli::cli_h1("Beginning clustering")
  }
  # Process the clustering queue
  while (length(clustering_queue) > 0) {
    fit <- clustering_queue[[1]]
    clustering_queue <- clustering_queue[-1]

    # If the collection contains only one network, add it to the final list
    if (fit$best_fit$M == 1) {
      list_model_binary <- append(list_model_binary, list(fit$best_fit))
      next
    }

    # Compute the dissimilarity matrix
    dist_bm <- compute_dissimilarity_matrix(collection = fit)
    # Partition the networks based on the dissimilarity matrix
    cl <- partition_networks_list_from_dissimilarity(
      networks_list = fit$A,
      dissimilarity_matrix = dist_bm,
      nb_groups = 2L
    )

    if (verbose) {
      cli::cli_h2("Trying to split the collection of {.val {fit$net_id}}")
    }
    # Fit models for the sub-collections
    fits <- future.apply::future_lapply(
      c(1, 2),
      function(k) {
        Z_init <- lapply(
          seq_along(fit$model_list),
          function(q) {
            if (!is.null(fit$model_list[[q]])) {
              return(fit$model_list[[q]]$Z[cl == k])
            } else {
              return(NULL)
            }
          }
        )
        dim(Z_init) <- c(fit$global_opts$Q1_max, fit$global_opts$Q2_max)

        if (verbose) {
          cli::cli_alert_info("Fitting a sub collection with : {.val {fit$net_id[cl == k]}}")
        }

        filtered_sep_BiSBM <- vector("list")
        filtered_sep_BiSBM$model <- fit$sep_BiSBM$model[cl == k]
        filtered_sep_BiSBM$BICL <- fit$sep_BiSBM$BICL[cl == k]
        filtered_sep_BiSBM$Z <- fit$sep_BiSBM$Z[cl == k]

        return(
          estimate_colBiSBM(
            netlist = fit$A[cl == k],
            colsbm_model = colsbm_model,
            net_id = fit$net_id[cl == k],
            distribution = distribution,
            nb_run = min(sum(cl == k), nb_run),
            Z_init = Z_init,
            global_opts = global_opts,
            fit_opts = fit_opts,
            sep_BiSBM = filtered_sep_BiSBM
          )
        )
      },
      future.seed = TRUE
    )

    # Decide whether to continue splitting or add to final list
    if (full_inference || (fits[[1]]$best_fit$BICL + fits[[2]]$best_fit$BICL > fit$best_fit$BICL)) {
      clustering_queue <- append(clustering_queue, fits)
      if (verbose) {
        cli::cli_alert_success("Splitting collections improved the BIC-L criterion")
      }
    } else {
      list_model_binary <- append(list_model_binary, list(fit$best_fit))
      if (verbose) {
        cli::cli_alert_danger("Splitting collections {.emph decreased} the BIC-L criterion")
      }
    }
  }

  # Final message indicating the end of clustering
  if (verbose) {
    cli::cli_alert_success("Finished clustering")
  }
  invisible(list_model_binary)
}

## Implement all of this in one big function (clusterize_networks)
## That can be used for both types of colSBM
# Reworked clustering procedure
# 1. Fit the full collection (if separated preferred suggest a cut ?)
# 2. Compute dissimilarity matrix
# 3. Perform a cut based on diss matrix (let user choose the number of sub-partitions ?)
# 4. Compute the sum of BICL and compare to the full collection
# If the sum of BICL > BICL:
# 5.a. Go back to step 2 for all sub collections
# 5.b This partition is currently the best


#' Computes a dissimilarity matrix between networks of the collections
#'
#' @param collection A bmpop or bisbmpop object on which to build the
#' dissimilarity matrix
#'
#' @param weight The weighting to apply to the block proportions. One of "max"
#' or "mean", defaults to "max".
#'
#' @param norm The norm to use, either one of "L1" or "L2". Defaults to "L2".
#'
#'  @return A matrix of size \eqn{M * M} containing the dissimilarity matrix
#' between the networks.
#'
#' @export
compute_dissimilarity_matrix <- function(
    collection,
    weight = "max",
    norm = "L2") {
  stopifnot(
    "Can't build the distance matrix, this is not a bmpop or bisbmpop object" =
      (inherits(collection, "bisbmpop") | inherits(collection, "bmpop"))
  )

  if (inherits(collection, "bisbmpop")) {
    dist_matrix <- compute_dissimilarity_matrix.bisbmpop(collection,
      weight = weight,
      norm = norm
    )
  }
  if (inherits(collection, "bmpop")) {
    dist_matrix <- compute_dissimilarity_matrix.bmpop(collection,
      weight = weight,
      norm = norm
    )
  }

  return(dist_matrix)
}

#' @export
compute_dissimilarity_matrix.bisbmpop <- function(
    collection,
    weight = "max",
    norm = "L2") {
  M <- collection$M
  dist_matrix <- matrix(0, nrow = M, ncol = M)

  dist_matrix <- outer(
    1:M, 1:M,
    Vectorize(function(i, j) {
      if (i == j) {
        return(0)
      }
      pis <- lapply(collection$best_fit$pim[c(i, j)], function(list) list[[1]])
      rhos <- lapply(collection$best_fit$pim[c(i, j)], function(list) list[[2]])

      dist_bisbmpop_max(
        pi = pis, rho = rhos,
        alpha = collection$best_fit$alpham[c(i, j)],
        weight = weight,
        norm = norm
      )
    })
  )

  return(dist_matrix)
}

#' @export
compute_dissimilarity_matrix.bmpop <- function(
    collection,
    weight = "max",
    norm = "L2") {
  M <- collection$M
  dist_matrix <- matrix(0, nrow = M, ncol = M)

  dist_matrix <- outer(
    1:M, 1:M,
    Vectorize(function(i, j) {
      if (i == j) {
        return(0)
      }
      pis <- lapply(collection$best_fit$pim[c(i, j)], function(list) list)

      dist_bmpop_max(
        pi = pis,
        alpha = collection$best_fit$alpham[c(i, j)],
        delta = collection$best_fit$delta[c(i, j)],
        directed = collection$directed,
        weight = weight,
        norm = norm
      )
    })
  )


  return(dist_matrix)
}

#' Partition networks according to the dissimilarity matrix computed with
#' `compute_dissimilarity_matrix`.
#'
#' @details This functions partition a provided networks list and outputs back
#' a vector
#'
#' @param networks_list The list of networks to partition
#' @param nb_groups An integer, the number of groups. Defaults to 2
#'
#' @return A vector, eventually named according to `networks_list` names.
#'
#' @export
partition_networks_list_from_dissimilarity <- function(
    networks_list,
    dissimilarity_matrix,
    nb_groups = 2L) {
  # Sanity checks
  check_networks_list(networks_list)
  check_dissimilarity_matrix(dissimilarity_matrix)
  M <- length(networks_list)

  # Checking correspondance of args
  if (nrow(dissimilarity_matrix) != M || ncol(dissimilarity_matrix) != M) {
    cli::cli_abort(c("{.arg dissimilarity_matrix} has incorrect dimensions.",
      "i" = "It should be ({toString(c(M,M))}) to match {.arg networks_list} length.",
      "x" = "And it is of size ({toString(dim(dissimilarity_matrix))})"
    ), call = rlang::caller_env())
  }


  if (M >= 3L) {
    # If there is more than 3 networks they are splitted using partition around
    # K-medioids
    cl <- cluster::pam(
      x = sqrt(dissimilarity_matrix),
      k = nb_groups,
      diss = TRUE,
      cluster.only = TRUE
    )
  } else {
    cl <- c(1L, 2L)
  }

  names(cl) <- names(networks_list)
  return(cl)
}

#' Perform the clustering procedure
#'
#' @param networks_list A list of matrices representing the networks to
#' partition
#' @inheritParams estimate_colBiSBM
#' @param full_inference A boolean indicating if we should do the procedure
#' until all networks are in their own partition.
#' @importFrom fs file_temp
#' @param save_path A file path specifying where to save the current clustering.
#' Useful for long running clustering. Defaults to
#' `fs::file_temp(pattern = "clustering", ext = "Rds")`.
#'
#' @export
clusterize_bipartite_collection <- function(networks_list,
                                            colsbm_model,
                                            net_id = NULL,
                                            distribution = "bernoulli",
                                            nb_run = 3L,
                                            global_opts = list(),
                                            fit_opts = list(),
                                            fit_init = NULL,
                                            full_inference = FALSE,
                                            save_path = fs::file_temp(
                                              pattern = "clustering",
                                              ext = "Rds"
                                            )) {
  # Sanity checks
  check_networks_list(networks_list = networks_list)
  check_bipartite_colsbm_models(colsbm_model = colsbm_model)
  check_is_integer_over_thresh(nb_run, thresh = 1L)
  check_colsbm_emission_distribution(emission_distribution = distribution)
  check_networks_list_match_emission_distribution(networks_list = networks_list, emission_distribution = distribution)

  global_opts <- utils::modifyList(default_global_opts_bipartite(netlist = networks_list), global_opts)

  # Starting the clustering
  cli::cli_h1(text = "Clustering on {.emph M = {length(networks_list)}} networks")

  # Initialization of clustering history list
  step <- 1L
  max_steps <- ceiling(log2(length(networks_list))) + 1L
  clustering_history_list <- list()
  # Initialization of needed values
  current_bicl <- NULL
  current_partition_bicl <- -Inf
  current_partition_indices <- rep(1, length(networks_list))

  has_bicl_increased_with_partition <- TRUE
  can_we_split_further <- TRUE

  while (full_inference ||
    (has_bicl_increased_with_partition && any(can_we_split_further))) {
    current_bicl <- append(current_bicl, current_partition_bicl)
    current_partition_matrices <- lapply(
      unique(current_partition_indices),
      function(k) {
        networks_list[current_partition_indices == k]
      }
    )
    clustering_history_list <- append(clustering_history_list, list(current_partition_matrices))
    # Beginning the fit of the full collection
    cli::cli_h2("Step {step}/{max_steps}")

    # max_substep <- length(current_partition_matrices)

    current_partition_fit <- lapply(seq_along(current_partition_matrices), function(idx) {
      current_substep <- idx
      current_M <- length(current_partition_matrices[[idx]])
      cli::cli_h3(
        "Sub-step {current_substep}/{length(current_partition_matrices)} - Fitting colSBM on {current_M} network{?s}."
      )

      estimate_colBiSBM(
        netlist = current_partition_matrices[[idx]],
        colsbm_model = colsbm_model,
        distribution = distribution,
        nb_run = nb_run,
        global_opts = global_opts
      )
    })
    # Should the procedure continue
    #  I will need to keep vector of bicl to be able to determine if a
    # partition must be splitted or not
    current_partition_bicl <- sum(sapply(current_partition_fit, function(fit) fit$best_fit$BICL))
    has_bicl_increased_with_partition <- tail(current_bicl, 1) < current_partition_bicl
    can_we_split_further <- sapply(current_partition_fit, function(fit) fit$M) > 1L

    current_diss_matrix <- lapply(current_partition_fit, function(fit) compute_dissimilarity_matrix(collection = fit))
    current_partition_indices <- unlist(lapply(seq_along(current_partition_fit), function(fit_idx) {
      if (can_we_split_further[[fit_idx]]) {
        2 * (fit_idx - 1) + partition_networks_list_from_dissimilarity(
          networks_list = current_partition_fit[[fit_idx]]$A,
          dissimilarity_matrix = current_diss_matrix[[fit_idx]],
          nb_groups = 2L
        )
      } else {
        setNames(2 * (fit_idx - 1) + 1, current_partition_fit[[fit_idx]]$net_id)
      }
    }))
    step <- step + 1
  }
  return(current_partition_fit)
}

#' Convert to tree
#'
#' @importFrom phylogram read.dendrogram
#' @importFrom stringr str_replace_all
#'
#' @param clustering A nested list given by one of the clusterize function from
#' which to extract the clustering tree.
#' @param invalid_char_to_replace_regex A regex string used by
#' `stringr::str_replace_all()` to clean net_ids before they are processed.
#' @param net_id_width An integer to truncate long net_id and prevent messy
#' plots. Defaults to 20.
#'
#' @return A dendrogram object
#'
#' @export
#'
#' @details
#' This function converts the nested list given by the clusterize functions
#' in Newick tree format that is read by `phylogram::read.dendrogram`.
#'
#' The code is adapted from this StackOverflow answer :
#' https://stackoverflow.com/questions/45091691/convert-a-nested-list-in-dendrogram-tree-with-r
extract_clustering_dendrogram <- function(
    clustering,
    invalid_char_to_replace_regex = "\\(|\\)",
    net_id_width = 20L) {
  partition <- extract_best_partition(clustering, unnest = FALSE)
  if ((
    inherits(partition, "fitBipartiteSBMPop") ||
      inherits(partition, "fitSimpleSBMPop")) &&
    !is.list(partition)) {
    partition <- list(partition)
  }
  net_id_clustering <- rapply(partition,
    function(collection) {
      collection[["net_id"]]
    },
    how = "list"
  )

  net_id_clustering <- rapply(net_id_clustering, function(vec_net_id) {
    vec_net_id <- stringr::str_replace_all(vec_net_id, invalid_char_to_replace_regex, "")
    vec_net_id <- stringr::str_trunc(string = vec_net_id, width = net_id_width)
  }, how = "list")

  newick_str <- paste0(lapply(net_id_clustering, function(y) paste0("(", paste0(y, collapse = ","), ")")), collapse = ",")

  # remove unwanted characters
  newick_str <- gsub('\"|c|list| ', "", newick_str)
  newick_str <- paste0("(", newick_str, ");")

  # remove brackets from single term list object
  newick_str <- stringr::str_replace_all(newick_str, "\\([a-z]*\\)", function(x) gsub("^\\(|\\)$", "", x))

  return(phylogram::read.dendrogram(text = newick_str))
}

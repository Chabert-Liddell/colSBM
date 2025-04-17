#' Partition of a collection of unipartite networks based on their common
#' mesoscale structures
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "delta", "deltapi",
#' @param directed A boolean, should the networks be considered as directed or not.
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
#' @param temp_save_path A string, the path where to save the temporary results.
#' Defaults to a temporary file.
#'
#' @importFrom future.apply future_lapply
#' @import cli
#' @importFrom utils modifyList
#'
#' @return A list with two elements:
#' \item{partition}{A list of models giving the best partition.}
#' \item{cluster}{A vector of integers giving the cluster of each network.}
#'
#' @details
#' This functions make call to `estimate_colSBM`.
#' @export
#'
#' @seealso [colSBM::estimate_colSBM()],
#' \code{\link[colSBM]{fitSimpleSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#' alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
#' alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
#' first_collection <- generate_unipartite_collection(
#'   n = 50,
#'   pi = c(0.5, 0.5),
#'   alpha = alpha1, M = 2
#' )
#' second_collection <- generate_unipartite_collection(
#'   n = 50,
#'   pi = c(0.5, 0.5),
#'   alpha = alpha2, M = 2
#' )
#'
#' netlist <- append(first_collection, second_collection)
#'
#' \dontrun{
#' cl_separated <- clusterize_networks(
#'   netlist = netlist,
#'   colsbm_model = "iid"
#' )
#' }
clusterize_unipartite_networks <- function(netlist,
                                           colsbm_model,
                                           directed = FALSE,
                                           net_id = NULL,
                                           distribution = "bernoulli",
                                           nb_run = 3L,
                                           global_opts = list(),
                                           fit_opts = list(),
                                           fit_init = NULL,
                                           full_inference = FALSE,
                                           verbose = TRUE,
                                           temp_save_path = tempfile(fileext = ".Rds")) {
  check_unipartite_colsbm_models(colsbm_model = colsbm_model)
  # Check if a netlist is provided, try to cast it if not
  check_networks_list(networks_list = netlist)
  net_id <- check_net_id_and_initialize(net_id = net_id, networks_list = netlist)
  check_colsbm_emission_distribution(emission_distribution = distribution)
  check_networks_list_match_emission_distribution(
    networks_list = netlist,
    emission_distribution = distribution
  )
  go <- default_global_opts_unipartite(netlist = netlist)
  go$plot_details <- 0
  go <- utils::modifyList(go, global_opts)
  global_opts <- go
  fo <- default_fit_opts_unipartite()
  fo <- utils::modifyList(fo, fit_opts)
  fit_opts <- fo

  # Fit the initial model on the full collection
  if (verbose) {
    if (!is.null(temp_save_path)) {
      cli::cli_alert_info("A save file will be created at {.val {temp_save_path}} and updated after each step")
    }
    cli::cli_h1("Fitting the full collection")
  }
  start_time <- Sys.time()
  my_sbmpop <- estimate_colSBM(
    netlist = netlist,
    colsbm_model = colsbm_model,
    directed = directed,
    net_id = net_id,
    distribution = distribution,
    nb_run = nb_run,
    global_opts = global_opts,
    fit_opts = fit_opts
  )

  clustering_queue <- list(my_sbmpop)
  list_model_binary <- list()
  cluster <- rep(1, length(netlist))
  names(cluster) <- net_id

  cluster_history <- as.data.frame(matrix(cluster, nrow = 1L))



  if (verbose) {
    cli::cli_h1("Beginning clustering")
  }
  # Process the clustering queue
  while (length(clustering_queue) > 0) {
    if (!is.null(temp_save_path)) {
      saveRDS(list(
        clustering_queue = clustering_queue,
        list_model_binary = list_model_binary,
        cluster_history = cluster_history
      ), temp_save_path)
    }
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
          seq_along(fit$model_list[[1]]),
          function(q) {
            lapply(
              seq_along(fit$model_list[[1]][[q]]),
              function(j) fit$model_list[[1]][[q]][[j]]$Z[cl == k]
            )
          }
        )

        if (verbose) {
          cli::cli_alert_info("Fitting a sub collection with : {.val {fit$net_id[cl == k]}}")
        }

        return(
          estimate_colSBM(
            netlist = fit$A[cl == k],
            colsbm_model = colsbm_model,
            net_id = fit$net_id[cl == k],
            distribution = distribution,
            nb_run = min(sum(cl == k), nb_run),
            Z_init = Z_init,
            global_opts = global_opts,
            fit_opts = fit_opts,
            fit_sbm = fit$fit_sbm[cl == k],
          )
        )
      },
      future.seed = TRUE
    )


    bicl_increased <- (fits[[1]]$best_fit$BICL + fits[[2]]$best_fit$BICL > fit$best_fit$BICL)
    # Decide whether to continue splitting or add to final list
    if (full_inference || bicl_increased) {
      clustering_queue <- append(clustering_queue, fits)
      if (verbose && bicl_increased) {
        cli::cli_alert_success("Splitting collections improved the BIC-L criterion")
      }
      if (verbose && full_inference) {
        cli::cli_alert_info("Full inference mode enabled, continuing to split collections")
      }
      prev_cluster <- unique(cluster[fit$net_id])

      #  Making room for a new cluster
      cluster[cluster > prev_cluster] <- cluster[cluster > prev_cluster] + 1

      # Assign the new cluster to the networks
      cluster[fit$net_id[cl == 2]] <- prev_cluster + 1
      cluster_history <- rbind(cluster_history, matrix(unname(cluster), nrow = 1))
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

  colnames(cluster_history) <- net_id

  output_list <- list(
    partition = list_model_binary,
    cluster = cluster,
    elapsed_time = Sys.time() - start_time,
    cluster_history = cluster_history
  )
  if (!is.null(temp_save_path)) {
    saveRDS(output_list, temp_save_path)
    if (verbose) {
      cli::cli_alert_info("The final results are saved at {.val {temp_save_path}}")
    }
  }
  return(output_list)
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
#' @param temp_save_path A string, the path where to save the temporary results.
#' Defaults to a temporary file.
#'
#' @importFrom future.apply future_lapply
#' @import cli
#' @importFrom utils modifyList
#'
#' @return A list of models for the recursive partition of
#' the collection of networks.
#'
#' This functions make call to `estimate_colBiSBM`.
#' @export
#'
#' @seealso [colSBM::clusterize_unipartite_networks()], [colSBM::estimate_colBiSBM()],
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
                                          verbose = TRUE,
                                          temp_save_path = tempfile(fileext = ".Rds")) {
  check_bipartite_colsbm_models(colsbm_model = colsbm_model)
  # Check if a netlist is provided, try to cast it if not
  check_networks_list(networks_list = netlist)
  net_id <- check_net_id_and_initialize(net_id = net_id, networks_list = netlist)
  check_colsbm_emission_distribution(emission_distribution = distribution)
  check_networks_list_match_emission_distribution(
    networks_list = netlist,
    emission_distribution = distribution
  )
  go <- default_global_opts_bipartite(netlist = netlist)
  go <- utils::modifyList(go, global_opts)
  global_opts <- go
  fo <- default_fit_opts_bipartite()
  fo <- utils::modifyList(fo, fit_opts)
  fit_opts <- fo

  # Fit the initial model on the full collection
  if (verbose) {
    if (!is.null(temp_save_path)) {
      cli::cli_alert_info("A save file will be created at {.val {temp_save_path}} and updated after each step")
    }
    cli::cli_h1("Fitting the full collection")
  }
  start_time <- Sys.time()
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
  cluster <- rep(1, length(netlist))
  names(cluster) <- net_id
  clustering_history <- as.data.frame(matrix(cluster, nrow = 1L))

  if (verbose) {
    cli::cli_h1("Beginning clustering")
  }
  # Process the clustering queue
  while (length(clustering_queue) > 0) {
    if (!is.null(temp_save_path)) {
      saveRDS(list(
        clustering_queue = clustering_queue,
        list_model_binary = list_model_binary,
        clustering_history = clustering_history
      ), temp_save_path)
    }
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


    bicl_increased <- (fits[[1]]$best_fit$BICL + fits[[2]]$best_fit$BICL > fit$best_fit$BICL)
    # Decide whether to continue splitting or add to final list
    if (full_inference || bicl_increased) {
      clustering_queue <- append(clustering_queue, fits)
      if (verbose && bicl_increased) {
        cli::cli_alert_success("Splitting collections improved the BIC-L criterion")
      }
      if (verbose && full_inference) {
        cli::cli_alert_info("Full inference mode enabled, continuing to split collections")
      }
      prev_cluster <- unique(cluster[fit$net_id])

      #  Making room for a new cluster
      cluster[cluster > prev_cluster] <- cluster[cluster > prev_cluster] + 1

      # Assign the new cluster to the networks
      cluster[fit$net_id[cl == 2]] <- prev_cluster + 1
      clustering_history <- rbind(clustering_history, matrix(unname(cluster), nrow = 1))
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
  colnames(clustering_history) <- net_id
  output_list <- list(
    partition = list_model_binary,
    cluster = cluster,
    elapsed_time = Sys.time() - start_time,
    clustering_history = clustering_history
  )
  if (!is.null(temp_save_path)) {
    saveRDS(output_list, temp_save_path)
    if (verbose) {
      cli::cli_alert_info("The final results are saved at {.val {temp_save_path}}")
    }
  }
  return(output_list)
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
#' @return A matrix of size \eqn{M * M} containing the dissimilarity matrix
#' between the networks.
#'
#' @keywords internal
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

#' Compute the dissimilarity matrix for a collection of bipartite networks
#'
#' @inheritParams compute_dissimilarity_matrix
#' @keywords internal
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

#' Compute the dissimilarity matrix for a collection of networks
#' @inheritParams compute_dissimilarity_matrix
#' @keywords internal
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
#' @param dissimilarity_matrix The dissimilarity matrix computed with
#' `compute_dissimilarity_matrix`.
#' @param nb_groups An integer, the number of groups. Defaults to 2
#'
#' @return A vector, eventually named according to `networks_list` names.
#'
#' @keywords internal
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
    cl <- cutree(hclust(as.dist(sqrt(dissimilarity_matrix)), method = "single"), k = nb_groups)
  } else {
    cl <- c(1L, 2L)
  }

  names(cl) <- names(networks_list)
  return(cl)
}

#' Clusterize bipartite networks in a bottom up fashion with graphon distance
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colSBM to use, one of "iid", "pi", "rho", "pirho",
#' "delta", "deltapi".
#' @param net_id A vector of string, the name of the networks.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param nb_run An integer, the number of run the algorithm do. Defaults to 3.
#' @param global_opts Global options for the outer algorithm and the output
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init WIP A list of fitted collections from which to start the
#' fusions
#' Optional fit init from where initializing the algorithm.
#' @param full_inference The default "FALSE", the algorithm stop once splitting
#' groups of networks does not improve the BICL criterion. If "TRUE", then
#' continue to split groups until a trivial classification of one network per
#' group.
#' @param keep_history A boolean, should the function keep the history of the
#' fusions or not. Default to FALSE to reduce file size.
#' Note that the fusion history is saved in the temporary file.
#' @param verbose A boolean, should the function be verbose or not. Default to
#' TRUE.
#' @param temp_save_path A string, the path where to save the temporary results.
#' Defaults to a temporary file.
#'
#' @export
#'
#' @import cli
#' @import progressr
#' @importFrom stats as.dist
#'
#' @return A list of models for the recursive partition of
#' the collection of networks.
clusterize_bipartite_networks_graphon <- function(
    netlist,
    colsbm_model,
    net_id = NULL,
    distribution = "bernoulli",
    nb_run = 3L,
    global_opts = list(),
    fit_opts = list(),
    fit_init = NULL, # Use this to store a list of fits from which to start clustering
    full_inference = FALSE,
    keep_history = FALSE,
    verbose = TRUE,
    temp_save_path = tempfile(fileext = ".Rds")) {
  # Check the colSBM model
  check_bipartite_colsbm_models(colsbm_model = colsbm_model)
  stopifnot("Only iid model is implemented" = colsbm_model == "iid")
  # Check if a netlist is provided, try to cast it if not
  check_networks_list(networks_list = netlist)
  net_id <- check_net_id_and_initialize(net_id = net_id, networks_list = netlist)
  check_colsbm_emission_distribution(emission_distribution = distribution)
  check_networks_list_match_emission_distribution(
    networks_list = netlist,
    emission_distribution = distribution
  )
  go <- default_global_opts_bipartite(netlist = netlist)
  go <- utils::modifyList(go, global_opts)
  global_opts <- go
  fo <- default_fit_opts_bipartite()
  fo <- utils::modifyList(fo, fit_opts)
  fit_opts <- fo

  if (verbose) {
    cli::cli_text("Temporary results will be saved as {.val {temp_save_path}} after each step")
  }
  if (is.null(fit_init)) {
    # Initializing separated collections
    cli::cli_h1("Beginning the clustering")
    cli::cli_h2("Fitting separated BiSBM models")
    p <- progressr::progressor(along = netlist)
    collections <- colsbm_lapply(seq_along(netlist), function(i) {
      col <- estimate_colBiSBM(
        netlist = list(netlist[[i]]),
        net_id = c(net_id[[i]]),
        colsbm_model = colsbm_model,
        nb_run = nb_run,
        distribution = distribution,
        global_opts = global_opts,
        fit_opts = fit_opts
      )
      if (exists("p")) {
        p(sprintf("Fitted network %s", paste0(col$net_id)))
      }
      col
    },
    backend = global_opts$backend,
    nb_cores = global_opts$nb_cores
    )
  } else {
    collections <- fit_init
  }

  cluster <- seq_len(length(netlist))
  names(cluster) <- net_id

  if (!is.null(fit_init) & verbose) {
    cli::cli_alert_info("Starting from a list of fits")
    cli::cli_alert_warning("This feature is still experimental and cluster vector may not be accurate\n\n")
  }

  compute_bicl_partition <- function(collections) {
    return(sum(sapply(collections, function(col) col$best_fit$BICL)))
  }

  # Historique des fusions
  if (verbose) {
    cli::cli_alert_info("Fully separated BIC-L is {.val {compute_bicl_partition(collections)}}")
  }
  bicl_history <- c(compute_bicl_partition(collections))
  fusion_history <- list(collections)

  # Fonction pour calculer les distances de graphon entre toutes les paires de collections
  compute_distances <- function(collections) {
    M <- length(collections)
    dist_matrix <- matrix(Inf, nrow = M, ncol = M)
    for (i in seq_len(M)) {
      for (j in seq_len(i - 1)) {
        dist_matrix[i, j] <- dist_graphon_bipartite_all_permutations(
          pis = list(collections[[i]]$best_fit$parameters$pi[[1]], collections[[j]]$best_fit$parameters$pi[[1]]),
          rhos = list(collections[[i]]$best_fit$parameters$rho[[1]], collections[[j]]$best_fit$parameters$rho[[1]]),
          alphas = list(collections[[i]]$best_fit$parameters$alpha, collections[[j]]$best_fit$parameters$alpha)
        )
      }
    }
    return(dist_matrix)
  }

  # Function to generate and sort pairs of indices based on distance matrix
  generate_sorted_pairs <- function(dist_matrix) {
    pairs <- which(lower.tri(dist_matrix), arr.ind = TRUE)
    distances <- dist_matrix[lower.tri(dist_matrix)]
    sorted_indices <- order(distances)
    sorted_pairs <- pairs[sorted_indices, ]
    return(sorted_pairs)
  }

  max_steps <- length(collections) - 1
  step <- 1
  has_bicl_increased <- TRUE
  # Loop to merge collections
  while (length(collections) > 1 && (has_bicl_increased || full_inference)) {
    saveRDS(list(fusion_history = fusion_history, bicl_history = bicl_history), temp_save_path)
    if (verbose) {
      cli::cli_h2("Step {.val {step}} on a max of {.val {max_steps}} step{?s}")
      cli::cli_alert_info("Computing distances between collections")
    }
    dist_matrix <- compute_distances(collections)
    sorted_pairs <- generate_sorted_pairs(dist_matrix)
    sorted_pairs <- matrix(sorted_pairs, ncol = 2L)
    i <- sorted_pairs[1, 1]
    j <- sorted_pairs[1, 2]
    if (verbose) {
      cli::cli_alert_info("Try merging{cli::qty({collections[[i]]$net_id})} network{?s} {.val {collections[[i]]$net_id}} with{cli::qty({collections[[j]]$net_id})} network{?s} {.val {collections[[j]]$net_id}}")
    }
    candidate_collection <- estimate_colBiSBM(
      netlist = c(collections[[i]]$A, collections[[j]]$A),
      net_id = c(collections[[i]]$net_id, collections[[j]]$net_id),
      colsbm_model = colsbm_model,
      nb_run = nb_run,
      distribution = distribution,
      global_opts = global_opts,
      fit_opts = fit_opts
    )

    new_collection <- candidate_collection
    remaining_collections <- collections[-c(i, j)]
    has_bicl_increased <-
      # Test if there are remaining collections and compute their BICL or output 0
      ifelse(length(remaining_collections) >= 1,
        compute_bicl_partition(collections[-c(i, j)]),
        0
      ) +
        new_collection$best_fit$BICL > bicl_history[step]
    if (has_bicl_increased || full_inference) {
      # Mettre à jour les collections
      collections <- collections[-c(i, j)]
      collections <- c(collections, list(new_collection))
      bicl_history <- c(bicl_history, compute_bicl_partition(collections))
      cli::cli_alert_info("Current collection BIC-L is {.val {compute_bicl_partition(collections)}}")
      if (!has_bicl_increased && full_inference && verbose) {
        cli::cli_alert_info("Full inference requested, clustering will continue, but BIC-L has not improved")
      }

      net_id_to_merge <- new_collection$net_id
      cluster_to_merge <- cluster[which(names(cluster) %in% net_id_to_merge)]

      cluster[which(names(cluster) %in% net_id_to_merge)] <- rep(min(cluster_to_merge), length(cluster_to_merge))

      cluster[cluster > max(cluster_to_merge)] <- cluster[cluster > max(cluster_to_merge)] - 1

      # Ajouter à l'historique des fusions
      fusion_history <- c(fusion_history, list(collections))
      step <- step + 1
      if (verbose) {
        cli::cli_alert_info("Collections updated, {length(collections)} remaining")
      }
    } else {
      if (verbose) {
        cli::cli_alert_info("No improvement in BIC-L, clustering will stop")
      }
    }
  }
  if (verbose) {
    cli::cli_h1("Clustering of bipartite networks completed")
  }
  out <- list(
    partition = tail(fusion_history, 1)[[1]],
    cluster = cluster,
    bicl_history = bicl_history
  )

  if (keep_history) {
    out[["fusion_history"]] <- fusion_history
  }

  return(out)
}

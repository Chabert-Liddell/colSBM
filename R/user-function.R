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
#' @param global_opts Global options for the outer algorithm and the output
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#'
#' @return A bmpop object listing a collection of models for the collection.
#' of networks
#' @export
#'
#' @seealso [colSBM::clusterize_networks()], \code{\link[colSBM]{bmpop}},
#' \code{\link[colSBM]{fitSimpleSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#' # Trivial example with Gnp networks:
#' Net <- lapply(list(.7, .7, .2, .2),
#'               function(p) {
#'                 A <- matrix(0, 15, 15 )
#'                A[lower.tri(A)][sample(15*14/2, size = round(p*15*14/2))] <- 1
#'                 A <- A + t(A)
#'               })
#' \dontrun{cl <- estimate_colSBM(Net,
#'                      colsbm_model = "delta",
#'                      directed = FALSE,
#'                      model = "bernoulli",
#'                      nb_run = 1
#'                      )}
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

    # go is used to temporarily store the default global_opts
    go <- list(Q_min = 1L,
                             Q_max = floor(log(sum(vapply(netlist, "nrow", .1))))+2,
                             sbm_init = TRUE,
                             spectral_init = TRUE,
                             nb_init = 10L,
                             nb_models = 5L,
                             depth = 1L,
                             plot_details = 1L,
                             max_pass = 10L,
                             verbosity = 0L,
                             nb_cores = 1L)
    go <- utils::modifyList(go, global_opts)
    global_opts <- go
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
      if(is.null(fit_sbm) && global_opts$sbm_init) {
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
            },
            mc.cores = nb_cores,
            mc.silent = TRUE
          )
      }
      # tmp_fits run nb_run times a full model selection procedure
      # (the one from the research paper)
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
          }, mc.progress = TRUE, mc.cores = min(nb_run,nb_cores),
          mc.stdout = "output"
        )
      my_bmpop <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$BICL,
                                             FUN.VALUE = .1))]]
      # Aggregate all runs for all model sizes
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

#' Estimate a colBiSBM on a collection of networks
#'
#' @param netlist A list of matrices.
#' @param colsbm_model Which colBiSBM to use, one of "iid", "pi", "rho",
#' "pirho", "deltapi".
#' @param net_id A vector of string, the name of the networks.
#' @param distribution A string, the emission distribution, either "bernoulli"
#' (the default) or "poisson"
#' @param nb_run An integer, the number of run the algorithm do.
#' @param global_opts Global options for the outer algorithm and the output.
#' See details.
#' @param fit_opts Fit options for the VEM algorithm
#' @param fit_init Do not use!
#' Optional fit init from where initializing the algorithm.
#' @param Z_init A bi-dimensional list of size Q1_max x Q2_max containing
#' for each value a list of two vectors of clusters memberships. Default to 
#' NULL.
#'
#' @details The list of parameters \code{global_opts} essentially tunes the
#' optimization process and the variational EM algorithm,
#' with the following parameters
#'  \itemize{
#'  \item{\code{nb_cores} }{integer for number of cores used for
#'  parallelization. Default is 1}
#'  \item{\code{verbosity} }{integer for verbosity (0, 1, 2, 3, 4). Default is 1.
#'   0 will disable completely the output of the function.}
#'  \item{\code{Q1_max} }{integer for the max size in row to explore. Default is
#'  computed with the following formula:
#' `floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2)`}
#'  \item{\code{Q2_max} }{integer for the max size in columns to explore. Default is
#'  computed with the following formula:
#' `floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2)`}
#'  \item{\code{nb_models} }{the number of models to keep for each values of Q1,Q2.
#'  Default is 5.}
#'  \item{\code{depth} }{specifies how large will the moving window be. Default is 1,
#'  meaning the window will go from (Q1 - 1, Q2 - 1) to (Q1 + 1, Q2 + 1) and all
#'  the values in the square defined.}
#'  \item{\code{plot_details} }{integer to control the display of the exploration and
#'  moving window process. Values are 0 or 1. Default is 1.}
#'  \item{\code{max_pass} }{the maximum number of moving window passes that will be
#'  executed. Default is 10.}
#' }
#'
#' @return A bisbmpop object listing a collection of models for the collection.
#' of networks
#' @export
#'
#' @seealso [colSBM::clusterize_networks()], \code{\link[colSBM]{bisbmpop}},
#' \code{\link[colSBM]{fitBipartiteSBMPop}}, `browseVignettes("colSBM")`
#'
#' @examples
#' # Trivial example with Gnp networks:
#' Net <- lapply(
#'   list(.7, .7, .2, .2),
#'   function(p) {
#'     A <- matrix(0, 15, 15)
#'     A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
#'     A <- A + t(A)
#'   }
#' )
#' \dontrun{
#' cl <- estimate_colSBM(Net,
#'   colsbm_model = "delta",
#'   distribution = "bernoulli",
#'   nb_run = 1
#' )
#' }
estimate_colBiSBM <-
  function(netlist,
           colsbm_model,
           net_id = NULL,
           distribution = "bernoulli",
           nb_run = 3L,
           global_opts = list(),
           fit_opts = list(),
           fit_init = NULL,
           Z_init = NULL) {
    switch(colsbm_model,
      "iid" = {
        free_density <- FALSE
        free_mixture_row <- FALSE
        free_mixture_col <- FALSE
      },
      "pi" = {
        free_density <- FALSE
        free_mixture_row <- TRUE
        free_mixture_col <- FALSE
      },
      "rho" = {
        free_density <- FALSE
        free_mixture_row <- FALSE
        free_mixture_col <- TRUE
      },
      "pirho" = {
        free_density <- FALSE
        free_mixture_row <- TRUE
        free_mixture_col <- TRUE
      },
      "delta" = {
        free_density <- TRUE
        free_mixture_row <- FALSE
        free_mixture_col <- FALSE
      },
      "deltapi" = {
        free_density <- TRUE
        free_mixture_row <- TRUE
        free_mixture_col <- FALSE
      },
      stop(
        "colsbm_model unknown.",
        " Must be one of iid, pi, rho, pirho, delta or deltapi"
      )
    )

    # go is used to temporarily store the default global_opts
    go <- list(
      Q1_min = 1L,
      Q2_min = 1L,
      Q1_max = floor(log(sum(sapply(netlist, function(A) nrow(A)))) + 2),
      Q2_max = floor(log(sum(sapply(netlist, function(A) ncol(A)))) + 2),
      nb_init = 10L,
      nb_models = 5L,
      depth = 1L,
      plot_details = 1L,
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

    # To warn the user about the verbosity and nb_cores
    if (global_opts$verbosity >= 3 && global_opts$nb_cores > 1) {
      cat(
        "\nDue to the parallelization, the message logs",
        "might be in confusing order.",
        " You may want to use only one core for readability"
      )
    } 

    if (!is.null(fit_init)) {
      bisbmpop <- fit_init
    } else {
      if (is.null(net_id)) {
        net_id <- seq_along(netlist)
      }
      # tmp_fits run nb_run times a full model selection procedure
      # (the one from the research paper)
      tmp_fits <-
        bettermc::mclapply(
          seq(nb_run),
          function(x) {
            # Computes the number of cores to allocate per run
            global_opts$nb_cores <- max(
              1L,
              floor(global_opts$nb_cores / nb_run)
            )
            tmp_fit <- bisbmpop$new(
              netlist = netlist,
              net_id = net_id,
              distribution = distribution,
              free_density = free_density,
              free_mixture_row = free_mixture_row,
              free_mixture_col = free_mixture_col,
              global_opts = global_opts,
              fit_opts = fit_opts,
              Z_init = Z_init
            )
            tmp_fit$optimize()
            return(tmp_fit)
          },
          mc.progress = TRUE,
          mc.cores = min(nb_run, nb_cores),
          mc.stdout = "output"
        )
      # We choose the the bisbmpop to receive the best_fit in sense of the BICL
      bisbmpop <- tmp_fits[[which.max(vapply(tmp_fits, function(fit) fit$best_fit$BICL,
        FUN.VALUE = .1
      ))]]

      # We perform the procedure only if nb_run > 1
      if (nb_run > 1) {
        # For each Q1 and Q2, we compare all
        for (q1 in global_opts$Q1_max) {
          for (q2 in global_opts$Q2_max) {
            if (!is.null(bisbmpop$model_list[[q1, q2]])) {
              # All the models for the current q1 and q2 are stored
              models_comparison <- # Note : this object is a list
                c(bisbmpop$model_list[[q1, q2]], lapply(
                  tmp_fits,
                  function(fit) fit$model_list[[q1, q2]]
                ))
              # The best in the sense of the BICL is chosen
              bisbmpop$model_list[[q1, q2]] <-
                models_comparison[which.max(
                  vapply(models_comparison, function(model) model$BICL)
                )]

              # The same procedure is applied for the
              discarded_models_comparison <-
                c(bisbmpop$discarded_model_list[[q1, q2]], unlist(lapply(
                  tmp_fits,
                  function(fit) fit$discarded_model_list[[q1, q2]]
                )))
              bisbmpop$discarded_model_list[[q1, q2]] <-
                discarded_models_comparison[which.max(
                  vapply(discarded_models_comparison, function(model) model$BICL)
                )]
            }
          }
        }

        # We now update the criteria and best fit
        bisbmpop$store_criteria_and_best_fit()
        # The discarded model list is truncated
        bisbmpop$truncate_discarded_model_list()
      }
      rm(tmp_fits)
      gc()
      # At the end we show the results
      if (global_opts$verbosity >= 1) {
        cat(
          "\nAfter merging the", nb_run, "model runs,",
          "the criteria are the following:\n"
        )
        bisbmpop$print_metrics()
      }
      bisbmpop$choose_joint_or_separated()
    }
    return(bisbmpop)
  }

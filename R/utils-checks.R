#' Check that the number is an integer over a certain threshold
#' @noRd
check_is_integer_over_thresh <- function(
    int,
    thresh,
    arg = rlang::caller_arg(int),
    call = rlang::caller_env()) {
  rlang::check_required(int, arg = arg, call = call)
  if (!rlang::is_integerish(int, n = 1, finite = TRUE)) {
    cli::cli_abort(c("{.arg {arg}} must be {.obj_type_friendly {1L}}.",
      "x" = "You've provided {.obj_type_friendly {int}}."
    ), call = call)
  }
  if (int < thresh) {
    cli::cli_abort(c("{.arg {arg}} must be at least {.val {thresh}}.",
      "x" = "You've provided {.val {int}}."
    ), call = call)
  }
  return(as.integer(int))
}

#' Check function for networks_list
#' @noRd
check_networks_list <- function(networks_list,
                                min_length = 2L,
                                arg = rlang::caller_arg(networks_list),
                                call = rlang::caller_env()) {
  rlang::check_required(networks_list,
    arg = arg,
    call = call
  )
  if (!rlang::is_list(networks_list) |
    rlang::is_empty(networks_list) |
    !all(sapply(networks_list, is.matrix))) {
    cli::cli_abort("{.arg {arg}} must be a list of matrices.",
      arg = arg,
      call = call
    )
  }
  if (length(networks_list) < min_length) {
    cli::cli_abort(
      c("{.arg {arg}} must be of length at least {min_length}.",
        "x" = "You've provided a list of length {.val {length(networks_list)}}."
      ),
      arg = arg,
      call = call
    )
  }
}

#' Check net_id and eventually intialize it
#' @noRd
#' @param net_id A character vector representing network IDs.
#' @param networks_list A list of networks.
#' @param arg The argument name for error messages (default is the name of `net_id`).
#' @param call The calling environment (default is the caller environment).
#'
#' @importFrom rlang is_named is_empty
#'
#' @return The network IDs
check_net_id_and_initialize <- function(net_id, networks_list, arg = rlang::caller_arg(net_id), call = rlang::caller_env()) {
  if (rlang::is_named(networks_list)) {
    net_id <- names(networks_list)
  }
  if (rlang::is_empty(net_id)) {
    net_id <- seq_along(networks_list)
  }
  check_net_id(net_id, networks_list, arg = arg, call = call)
  return(net_id)
}

#' Check function for dissimilarity matrix
#' @noRd
check_dissimilarity_matrix <- function(dissimilarity_matrix,
                                       arg = rlang::caller_arg(dissimilarity_matrix),
                                       call = rlang::caller_env()) {
  rlang::check_required(dissimilarity_matrix,
    arg = arg,
    call = call
  )
  if (!is.matrix(dissimilarity_matrix) || !is.numeric(dissimilarity_matrix)) {
    cli::cli_abort(
      c("{.arg {arg}} must be a dissimilarity matrix.",
        "x" = "You've supplied {.obj_type_friendly {dissimilarity_matrix}} where {.obj_type_friendly {matrix(1,2)}} was expected."
      ),
      arg = arg,
      call = call
    )
  }
}

#' Check function for bipartite colSBM models
#' @noRd
check_bipartite_colsbm_models <- function(
    colsbm_model,
    arg = rlang::caller_arg(colsbm_model),
    call = rlang::caller_env()) {
  rlang::check_required(colsbm_model,
    arg = arg,
    call = call
  )
  rlang::arg_match(
    arg = colsbm_model,
    values = c("iid", "pi", "rho", "pirho"),
    error_arg = arg,
    error_call = call
  )
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
  return(list(free_mixture_row = free_mixture_row, free_mixture_col = free_mixture_col))
}

#' Check function for bipartite colSBM models
#' @noRd
check_unipartite_colsbm_models <- function(
    colsbm_model,
    arg = rlang::caller_arg(colsbm_model),
    call = rlang::caller_env()) {
  rlang::check_required(colsbm_model,
    arg = arg,
    call = call
  )
  rlang::arg_match(
    arg = colsbm_model,
    values = c("iid", "pi", "delta", "deltapi"),
    error_arg = arg,
    error_call = call
  )
}

#' Check colSBM emission distribution
#' @noRd
#' @param emission_distribution A character string specifying the emission distribution
#' @param arg The name of the argument
#' @param call The environment where the function was called
#' @return The emission distribution
check_colsbm_emission_distribution <- function(
    emission_distribution,
    arg = rlang::caller_arg(emission_distribution),
    call = rlang::caller_env()) {
  rlang::check_required(emission_distribution,
    arg = arg,
    call = call
  )
  rlang::arg_match(
    arg = emission_distribution,
    values = c("poisson", "bernoulli"),
    error_arg = arg,
    error_call = call
  )
}

#' Check matrices match emission distribution
#' @noRd
#' @param networks_list A list of matrices
#' @param emission_distribution A character string specifying the emission distribution
#' @param arg The name of the argument
#' @param call The environment where the function was called
#' @return The list of matrices
check_networks_list_match_emission_distribution <- function(
    networks_list,
    emission_distribution,
    arg = rlang::caller_arg(networks_list),
    distrib_arg = rlang::caller_arg(emission_distribution),
    call = rlang::caller_env()) {
  check_colsbm_emission_distribution(emission_distribution,
    arg = distrib_arg,
    call = call
  )

  switch(emission_distribution,
    poisson = {
      if (!all(sapply(networks_list, rlang::is_integerish))) {
        cli::cli_abort(
          c("For Poisson emission distribution, all matrices in {.arg {arg}} must have non-negative integer entries.",
            "x" = "You've provided a {.obj_type_friendly {networks_list[[1]]}} with non integer entries."
          ),
          call = call
        )
      }
      if (!all(sapply(networks_list, function(x) all(x >= 0)))) {
        cli::cli_abort(
          c("For Poisson emission distribution, all matrices in {.arg {arg}} must have non-negative integer entries.",
            "x" = "You've provided a matrix with negative entries."
          ),
          call = call
        )
      }
    },
    bernoulli = {
      if (!all(sapply(networks_list, function(x) all(rlang::is_integerish(x))))) {
        cli::cli_abort(
          c("Non integer entries in {.arg {arg}} are not allowed for Bernoulli emission distribution.",
            "x" = "You've provided a matrix with non integer entries, {.obj_type_friendly {networks_list[[1]]}}."
          )
        )
      }
      if (!all(sapply(networks_list, function(x) all(x %in% c(0, 1, NA))))) {
        cli::cli_abort("All matrices in {.arg {arg}} must have entries that are either 0 or 1.",
          call = call
        )
      }
    }
  )
}

#' Check Net ID
#'
#' This function checks if the `net_id` is a character vector and if its length matches the length of `networks_list`.
#'
#' @param net_id A character vector representing network IDs.
#' @param networks_list A list of networks.
#' @param arg The argument name for error messages (default is the name of `net_id`).
#' @param call The calling environment (default is the caller environment).
#' @return Throws an error if the checks fail.
#'
#' @importFrom rlang is_character is_integerish are_na
check_net_id <- function(net_id, networks_list, arg = rlang::caller_arg(net_id), call = rlang::caller_env()) {
  if (!(rlang::is_character(net_id) ||
    rlang::is_integerish(net_id, finite = TRUE)) ||
    any(rlang::are_na(net_id))) {
    cli::cli_abort("{.arg {arg}} must be a character or an integer vector.",
      arg = arg,
      call = call
    )
  }
  if (length(net_id) != length(networks_list)) {
    cli::cli_abort(
      c("{.arg {arg}} must have the same length as {.arg networks_list}.",
        "i" = "You've provided a vector of length {.val {length(net_id)}} where length {.val {length(networks_list)}} was expected."
      ),
      arg = arg,
      call = call
    )
  }
  if (anyDuplicated(net_id) > 0L) {
    cli::cli_abort("{.arg {arg}} must have unique values.",
      arg = arg,
      call = call
    )
  }
  return(net_id)
}

#' Check Backend
#' @noRd
#' @param backend A character string specifying the backend
#' @param arg The name of the argument
#' @param call The environment where the function was called
#' @return The backend
check_backend <- function(
    backend,
    arg = rlang::caller_arg(backend),
    call = rlang::caller_env()) {
  rlang::check_required(backend,
    arg = arg,
    call = call
  )
  rlang::arg_match(
    arg = backend,
    values = c("no_mc", "parallel", "future"),
    error_arg = arg,
    error_call = call
  )
}

#' Check Global Options
#'
#' This function checks if `global_opts` is a list and if `nb_cores` (if provided) is an integer greater than a threshold.
#'
#' @param global_opts A list of global options.
#' @param arg The argument name for error messages (default is the name of `global_opts`).
#' @param call The calling environment (default is the caller environment).
#' @return Throws an error if the checks fail.
check_global_opts <- function(global_opts, arg = rlang::caller_arg(global_opts), call = rlang::caller_env()) {
  rlang::check_required(global_opts, arg = arg, call = call)
  if (!is.list(global_opts)) {
    cli::cli_abort("{.arg {arg}} must be a list.",
      arg = arg,
      call = call
    )
  }
  if (!is.null(global_opts$Q1_min)) {
    check_is_integer_over_thresh(global_opts$Q1_min, thresh = 1L)
  }
  if (!is.null(global_opts$Q1_max)) {
    check_is_integer_over_thresh(global_opts$Q1_max, thresh = 1L)
  }
  if (!is.null(global_opts$Q2_min)) {
    check_is_integer_over_thresh(global_opts$Q2_min, thresh = 1L)
  }
  if (!is.null(global_opts$Q2_max)) {
    check_is_integer_over_thresh(global_opts$Q2_max, thresh = 1L)
  }
  if (!is.null(global_opts$nb_init)) {
    check_is_integer_over_thresh(global_opts$nb_init, thresh = 1L)
  }
  if (!is.null(global_opts$nb_models)) {
    check_is_integer_over_thresh(global_opts$nb_models, thresh = 1L)
  }
  if (!is.null(global_opts$backend)) {
    check_backend(backend = global_opts$backend, arg = "global_opts$backend", call = call)
  }
  if (!is.null(global_opts$plot_details)) {
    check_is_integer_over_thresh(global_opts$plot_details, thresh = 0L)
  }
  if (!is.null(global_opts$depth)) {
    check_is_integer_over_thresh(global_opts$depth, thresh = 1L)
  }
  if (!is.null(global_opts$max_pass)) {
    check_is_integer_over_thresh(global_opts$max_pass, thresh = 1L)
  }
  if (!is.null(global_opts$verbosity)) {
    check_is_integer_over_thresh(global_opts$verbosity, thresh = 0L)
  }
  if (!is.null(global_opts$nb_cores)) {
    check_is_integer_over_thresh(global_opts$nb_cores, thresh = 1L)
  }
}

#' Check Fit Options
#'
#' This function checks if `fit_opts` is a list.
#'
#' @param fit_opts A list of fit options.
#' @param arg The argument name for error messages (default is the name of `fit_opts`).
#' @param call The calling environment (default is the caller environment).
#' @return Throws an error if the checks fail.
check_fit_opts <- function(fit_opts, arg = rlang::caller_arg(fit_opts), call = rlang::caller_env()) {
  rlang::check_required(fit_opts, arg = arg, call = call)
  if (!rlang::is_list(fit_opts)) {
    cli::cli_abort("{.arg {arg}} must be a list.",
      arg = arg,
      call = call
    )
  }
  rlang::arg_match0(
    arg = fit_opts$algo_ve,
    values = c("fp"),
    arg_nm = rlang::caller_arg(arg),
    error_call = call
  )
  if (!rlang::is_bool(fit_opts$minibatch)) {
    cli::cli_abort("{.arg fit_opts$minibatch} must be a boolean.",
      call = call
    )
  }
  if (!is.null(fit_opts$verbosity)) {
    check_is_integer_over_thresh(fit_opts$verbosity, thresh = 0L)
  }
  if (!rlang::is_double(fit_opts$tolerance, finite = TRUE)) {
    cli::cli_abort("{.arg fit_opts$tolerance} must be a double.",
      call = call
    )
  }
  if (!is.null(fit_opts$greedy_exploration_max_steps)) {
    check_is_integer_over_thresh(
      fit_opts$greedy_exploration_max_steps,
      thresh = 1L,
      call = call
    )
  }
  if (!is.null(fit_opts$greedy_exploration_max_steps_without_improvement)) {
    check_is_integer_over_thresh(
      fit_opts$greedy_exploration_max_steps_without_improvement,
      thresh = 1L,
      call = call
    )
  }
}

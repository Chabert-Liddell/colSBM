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
      if (!all(sapply(networks_list, rlang::is_integerish, finite = TRUE))) {
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
      if (!all(sapply(networks_list, function(x) all(x %in% c(0, 1))))) {
        cli::cli_abort("All matrices in {.arg {arg}} must have entries that are either 0 or 1.",
          call = call
        )
      }
    }
  )
}

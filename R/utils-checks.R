#' Check that the number is an integer over a certain threshold
#' @noRd
check_is_integer_over_thresh <- function(
    int,
    thresh,
    arg = rlang::caller_arg(int),
    call = rlang::caller_env()) {
    rlang::check_required(int, arg = arg, call = call)
    if (!rlang::is_integerish(int, n = 1, finite = TRUE) || int < thresh) {
        cli::cli_abort(c("{.arg {arg}} must be {.obj_type_friendly {1L}} that is at least {.val {thresh}}.",
            "x" = "You've provided {.val {int}} as {.obj_type_friendly {int}}."
        ), call = call)
    }
    return(as.integer(int))
}

#' Check function for networks_list
#' @noRd
check_networks_list <- function(networks_list,
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
    if (rlang::has_length(x = networks_list, n = 1L)) {
        cli::cli_abort("{.arg {arg}} must be of length at least 2, not 1.",
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
    if (!is.matrix(dissimilarity_matrix)) {
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

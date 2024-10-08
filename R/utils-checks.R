#' Check function for networks_list
check_networks_list <- function(networks_list,
                                arg = rlang::caller_arg(networks_list),
                                call = rlang::caller_env()) {
    rlang::check_required(networks_list, call = call)
    if (!rlang::is_list(networks_list) |
        rlang::is_empty(networks_list) |
        !all(sapply(networks_list, is.matrix))) {
        cli::cli_abort("{.arg {arg}} must be a list of matrices.", call = call)
    }
    if (rlang::has_length(x = networks_list, n = 1L)) {
        cli::cli_abort(
            c("{.arg {arg}} must be of length at least 2, not 1.", "x" = "{}"),
            call = call
        )
    }
}

#' Check function for dissimilarity matrix
check_dissimilarity_matrix <- function(dissimilarity_matrix,
                                       arg = rlang::caller_arg(dissimilarity_matrix),
                                       call = rlang::caller_env()) {
    rlang::check_required(dissimilarity_matrix, call = call)
    if (!is.matrix(dissimilarity_matrix)) {
        cli::cli_abort(
            c("{.arg {arg}} must be a dissimilarity matrix.",
                "x" = "You've supplied {.obj_type_friendly {dissimilarity_matrix}} where {.obj_type_friendly {matrix(1,2)}} was expected."
            ),
            call = call
        )
    }
}

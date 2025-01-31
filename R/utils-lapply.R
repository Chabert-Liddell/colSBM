#' Function to use parallel computing with user-specified or OS-specific
#' parallel backend
#'
#' @param X The data
#' @param FUN The function to use in parallel
#' @param backend One of c("future", "parallel")
#' @param nb_cores Number of parallel cores for parallel
#' @param ...
#'
#' @return A list on which the FUN was applied with specified backend and
#' parameters
#' @noRd
#' @noMd
#'
#' @examples
#' future::plan(future::multisession)
#' X <- seq(1, 5)
#' colsbm_lapply(X, backend = "future", nb_cores = 2L)
colsbm_lapply <- function(X, FUN, backend = "future",
                          nb_cores = 1L,
                          ...) {
  if (!(backend %in% c("future", "parallel", "no_mc"))) {
    stop("Invalid backend. Choose 'parallel', 'future' or 'no_mc'.")
  }
  if (backend == "future") {
    stopifnot(
      "The 'future.apply' package must be installed and configured with a plan outside this function." = requireNamespace("future.apply", quietly = TRUE)
    )
    result <- future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
    return(result)
  }

  if (backend == "parallel") {
    stopifnot(
      "The 'parallel' package must be installed to use this backend." = requireNamespace("parallel", quietly = TRUE)
    )
    if (toupper(Sys.info()["sysname"]) == "WINDOWS") {
      cl <- parallel::makeCluster(nb_cores)
      result <- parallel::parLapply(cl, X, FUN, ...)
      parallel::stopCluster(cl)
    } else {
      result <- parallel::mclapply(X, FUN,
        mc.cores = nb_cores,
        mc.preschedule = FALSE
      )
    }
    return(result)
  }

  if (backend == "no_mc") {
    result <- lapply(X, FUN)
    return(result)
  }
}

#' Function to use parallel computing with user-specified or OS-specific
#' parallel backend
#'
#' @param X The data
#' @param FUN The function to use in parallel
#' @param backend One of c("future", "parallel", "bettermc")
#' @param nb_cores Number of parallel cores for parallel or bettermc
#' @param ...
#'
#' @return A list on which the FUN was applied with specified backend and
#' parameters
#' @noRd
#' @noMd
#'
#' @examples
#' X <- seq(1, 5)
#' colsbm_lapply(X, backend = "parallel", nb_cores = 2L)
colsbm_lapply <- function(X, FUN, backend = "parallel",
                          nb_cores = parallelly::availableCores(omit = 1L),
                          ...) {
  if (!(backend %in% c("future", "parallel", "bettermc", "no_mc"))) {
    stop("Invalid backend. Choose 'parallel', 'future', 'bettermc' or 'no_mc'.")
  }
  if (backend == "future") {
    if (!isNamespaceLoaded("future.apply")) {
      stop("The 'future.apply' package must be loaded and configured with a plan
           outside this function.")
    }
    result <- future.apply::future_lapply(X, FUN, ..., future.seed = NULL)
    return(result)
  }

  if (backend == "parallel") {
    if (!isNamespaceLoaded("parallel")) {
      stop("The 'parallel' package must be loaded to use this backend.")
    }
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

  if (backend == "bettermc") {
    if (!isNamespaceLoaded("bettermc")) {
      stop("The 'bettermc' package must be loaded to use this backend.")
    }
    result <- bettermc::mclapply(X, FUN, mc.cores = nb_cores)
    return(result)
  }
  if (backend == "no_mc") {
    result <- lapply(X, FUN)
    return(result)
  }
}

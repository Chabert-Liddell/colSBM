#' Compute the dissimilarity between 2 mesoscale structures
#'
#' @param pi A list of two probability vectors
#' @param alpha A list of two connectivity matrices
#' @param delta A vector of 2 density parameters (optional)
#' @param weight One of "max" (default) or "mean". See details
#' @param norm "L1"or "L2" norm for the computation
#' @param directed Are the structure of the networks directed?
#'
#' @return The dissimilarity between two mesoscale structure.
#'
#' @details If weight is "max" then the weight of each block is computed as
#'  \code{pmax(pi[[1]], pi[[2]])}. If "mean", then we take the average. "max"
#'  penalize to a greater extent the difference in block proportion between
#'  structure.
#'
#' @export
#'
#' @examples
#' pi <- list(c(0.5,0.5), c(0.1,0.9))
#' alpha <- list(matrix(c(0.9,0.1,
#'                        0.1, 0.05), byrow = TRUE, nrow = 2),
#'                matrix(c(0.9,0.1,
#'                        0.1, 0.05), byrow = TRUE, nrow = 2))
#' dist_bmpop_max(pi, alpha)
dist_bmpop_max <- function(pi, alpha, delta = c(1, 1), weight = "max",
                           norm = "L2", directed) {
  if (missing(directed)) {
    directed <- isSymmetric.matrix(alpha[[1]]) & isSymmetric.matrix(alpha[[2]])
  }
  if (!directed) {
    alpha[[1]][upper.tri(alpha[[1]])] <- 0
    alpha[[2]][upper.tri(alpha[[2]])] <- 0
  }
  if (missing(pi)) {
    d <- switch(norm,
      "L1" = sum(abs(alpha[[1]] / delta[1] - alpha[[2]] / delta[2])),
      "L2" = sum((alpha[[1]] / delta[1] - alpha[[2]] / delta[2])**2)
    )
  } else {
    w <- switch(weight,
      "max" = pmax(pi[[1]], pi[[2]]),
      "mean" = .5 * (pi[[1]] + pi[[2]])
    )
    d <- switch(norm,
      "L1" = sum(w %*%
        (abs(alpha[[1]] / delta[1] - alpha[[2]] / delta[2])) %*%
        w),
      "L2" = sum(w %*%
        ((alpha[[1]] / delta[1] - alpha[[2]] / delta[2])**2) %*%
        w)
    )
  }
  d
}

#' Compute the dissimilarity between 2 mesoscale structures for bipartite 
#' SBM
#'
#' @param pi A list of two probability vectors (row)
#' @param rho A list of two probability vectors (columns)
#' @param alpha A list of two connectivity matrices
#' @param delta A vector of 2 density parameters (optional)
#' @param weight One of "max" (default) or "mean". See details
#' @param norm "L1"or "L2" norm for the computation
#'
#' @return The dissimilarity between two mesoscale structure.
#'
#' @details If weight is "max" then the weight of each block is computed as
#'  \code{pmax(pi[[1]], pi[[2]])}. If "mean", then we take the average. "max"
#'  penalize to a greater extent the difference in block proportion between
#'  structure.
#'
#' @export
#'
#' @examples
#' pi <- list(c(0.5,0.5), c(0.1,0.9))
#' rho <- list(c(0.1,0.9), c(0.5,0.5))
#' alpha <- list(matrix(c(0.9,0.1,
#'                        0.1, 0.05), byrow = TRUE, nrow = 2),
#'                matrix(c(0.9,0.1,
#'                        0.1, 0.05), byrow = TRUE, nrow = 2))
#' dist_bisbmpop_max(pi, rho, alpha)
dist_bisbmpop_max <- function(
    pi,
    rho, alpha, delta = c(1, 1), weight = "max", norm = "L2") {
  if (missing(pi) || missing(rho)) {
    distance <- switch(norm,
      "L1" = sum(abs(alpha[[1]] / delta[1] - alpha[[2]] / delta[2])),
      "L2" = sum((alpha[[1]] / delta[1] - alpha[[2]] / delta[2])**2)
    )
  } else {
    w_pi <- switch(weight,
      "max" = pmax(pi[[1]], pi[[2]]),
      "mean" = .5 * (pi[[1]] + pi[[2]])
    )
    w_rho <- switch(weight,
      "max" = pmax(rho[[1]], rho[[2]]),
      "mean" = .5 * (rho[[1]] + rho[[2]])
    )
    distance <- switch(norm,
      "L1" = sum(as.vector(w_pi) %*%
        (abs(alpha[[1]] / delta[1] - alpha[[2]] / delta[2])) %*%
        as.vector(w_rho)),
      "L2" = sum(as.vector(w_pi) %*%
        ((alpha[[1]] / delta[1] - alpha[[2]] / delta[2])**2) %*%
        as.vector(w_rho))
    )
  }
  distance
}

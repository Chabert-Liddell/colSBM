

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
dist_bmpop_max <- function(pi, alpha, delta = c(1,1), weight = "max",
                           norm = "L2", directed) {
  if (missing(directed)) {
    directed <- isSymmetric.matrix(alpha[[1]]) & isSymmetric.matrix(alpha[[2]])
  }
  if (! directed) {
    alpha[[1]][upper.tri(alpha[[1]])] <- 0
    alpha[[2]][upper.tri(alpha[[2]])] <- 0
  }
  if (missing(pi)) {
    d <- switch(
      norm,
      "L1" = sum(abs(alpha[[1]]/delta[1] - alpha[[2]]/delta[2])),
      "L2" = sum((alpha[[1]]/delta[1] - alpha[[2]]/delta[2])**2)
    )
  } else {
    w <- switch(
      weight,
      "max" = pmax(pi[[1]], pi[[2]]),
      "mean" = .5*(pi[[1]] + pi[[2]])
    )
    d <- switch(
      norm,
      "L1" = sum(w %*%
                   (abs(alpha[[1]]/delta[1] - alpha[[2]]/delta[2])) %*%
                   w),
      "L2" = sum(w %*%
                   ((alpha[[1]]/delta[1] - alpha[[2]]/delta[2])**2) %*%
                   w)
    )
  }
  d
}


# dist_bmpop_mean <- function(pi, alpha, delta = c(1,1)) {
#   alpha_minus <-
#     lapply(seq(2),
#            function(m) {
#              alpha_m <- 0*alpha[[1]]
#              for (q in seq(nrow(alpha_m))) {
#                for (r in seq(ncol(alpha_m))) {
#                  pi_q <- pi[[m]][-q]/sum(pi[[m]][-q])
#                  pi_r <- pi[[m]][-r]/sum(pi[[m]][-r])
#                  alpha_m[q,r] <- pi_q %*% alpha[[m]][-q, -r] %*% pi_r
#                }
#              }
#              alpha_m
#            })
#   dist_mean <-
#     sum(pmin(pi[[1]], pi[[2]]) %*% ((alpha[[1]]/delta[1] - alpha[[2]]/delta[2])**2) %*%
#         pmin(pi[[1]], pi[[2]]))
#   for(q in nrow(alpham[[1]])) {
#     for (r in ncol(alpham[[2]])) {
#       dist_mean <- dist_mean +
#         max(0, pi[[1]][q] - pi[[2]][q])*max(0, pi[[1]][r] - pi[[2]][r])*
#         (a)
#     }
#   }
#     pmax(0, pi[[1]] - pi[[2]]) %*% ((alpha[[1]]/delta - alpha_minus[[2]]/delta)**2) %*%
#     pmax(0, pi[[1]] - pi[[2]]) +
#     pmax(0, pi[[2]] - pi[[1]]) %*% ((alpha[[2]]/delta - alpha_minus[[1]]/delta)**2) %*%
#     pmax(0, pi[[2]] - pi[[1]])
#
# }
#

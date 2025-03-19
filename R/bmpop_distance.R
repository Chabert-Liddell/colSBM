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
#' pi <- list(c(0.5, 0.5), c(0.1, 0.9))
#' alpha <- list(
#'   matrix(c(
#'     0.9, 0.1,
#'     0.1, 0.05
#'   ), byrow = TRUE, nrow = 2),
#'   matrix(c(
#'     0.9, 0.1,
#'     0.1, 0.05
#'   ), byrow = TRUE, nrow = 2)
#' )
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
#' pi <- list(c(0.5, 0.5), c(0.1, 0.9))
#' rho <- list(c(0.1, 0.9), c(0.5, 0.5))
#' alpha <- list(
#'   matrix(c(
#'     0.9, 0.1,
#'     0.1, 0.05
#'   ), byrow = TRUE, nrow = 2),
#'   matrix(c(
#'     0.9, 0.1,
#'     0.1, 0.05
#'   ), byrow = TRUE, nrow = 2)
#' )
#' dist_bisbmpop_max(pi, rho, alpha)
dist_bisbmpop_max <- function(
    pi,
    rho, alpha, delta = c(1, 1), weight = "max", norm = "L2") {
  if (missing(pi) || missing(rho)) {
    distance <- switch(norm,
      "L1" = sum(abs(alpha[[1]] - alpha[[2]])),
      "L2" = sum((alpha[[1]] - alpha[[2]])**2)
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
        (abs(alpha[[1]] - alpha[[2]])) %*%
        as.vector(w_rho)),
      "L2" = sum(as.vector(w_pi) %*%
        ((alpha[[1]] - alpha[[2]])**2) %*%
        as.vector(w_rho))
    )
  }
  distance
}

#' Graphon distance for bipartite SBM
#' @param pis A list of two probability vectors (row)
#' @param rhos A list of two probability vectors (columns)
#' @param alphas A list of two connectivity matrices
#'
#' @importFrom utils head tail
#'
#' @return The graphon distance between two mesoscale structure.
#' @details The graphon distance is computed as the L2 norm between the
#' graphons of the two structures. Please note that this does not take into
#' account the possible permutation of the blocks.
#' @keywords internal
graphon_distance_bipartite <- function(pis, rhos, alphas) {
  # Extracting to pi1 and pi2
  pi1 <- as.vector(pis[[1]])
  pi2 <- as.vector(pis[[2]])

  # Extracting to rho1 and rho2
  rho1 <- rhos[[1]]
  rho2 <- rhos[[2]]

  # Extracting to alpha1 and alpha2
  alpha1 <- matrix(alphas[[1]], nrow = length(pi1), ncol = length(rho1))
  alpha2 <- matrix(alphas[[2]], nrow = length(pi2), ncol = length(rho2))
  dimnames(alpha1) <- list(
    paste0("row1.", seq_len(nrow(alpha1))),
    paste0("col1.", seq_len(ncol(alpha1)))
  )
  dimnames(alpha2) <- list(
    paste0("row2.", seq_len(nrow(alpha2))),
    paste0("col2.", seq_len(ncol(alpha2)))
  )

  # Cumsums of pi1 pi2 rho1 rho2
  pi1_cumsum <- c(0, cumsum(pi1))
  names(pi1_cumsum) <- paste0("row1.", seq_len(length(pi1_cumsum)) - 1)
  pi2_cumsum <- c(0, cumsum(pi2))
  names(pi2_cumsum) <- paste0("row2.", seq_len(length(pi2_cumsum)) - 1)
  rho1_cumsum <- c(0, cumsum(rho1))
  names(rho1_cumsum) <- paste0("col1.", seq_len(length(rho1_cumsum)) - 1)
  rho2_cumsum <- c(0, cumsum(rho2))
  names(rho2_cumsum) <- paste0("col2.", seq_len(length(rho2_cumsum)) - 1)


  min_pis <- outer(tail(pi1_cumsum, -1), tail(pi2_cumsum, -1), pmin)

  max_pis <- outer(head(pi1_cumsum, -1), head(pi2_cumsum, -1), pmax)

  # Here we compute the min(rho_l, rho_l')
  min_rhos <- outer(tail(rho1_cumsum, -1), tail(rho2_cumsum, -1), pmin)
  # Here we compute the max(rho_l-1, rho_l-1')
  max_rhos <- outer(head(rho1_cumsum, -1), head(rho2_cumsum, -1), pmax)

  # All possible surfaces
  inter_surf <- outer(pmax(min_pis - max_pis, 0), t(pmax(min_rhos - max_rhos, 0)))
  inter_surf <- aperm(inter_surf, c(1, 4, 2, 3)) # Adjusting dimensions

  # All possible connectivity diffrences
  diff_alpha <- outer(alpha1, alpha2, "-")

  # Compute the distance
  sum(inter_surf * (diff_alpha^2))
}

#' Graphon distance for bipartite SBM over all permutations of the blocks
#' @param pis A list of two probability vectors (row)
#' @param rhos A list of two probability vectors (columns)
#' @param alphas A list of two connectivity matrices
#'
#' @return The graphon distance between two mesoscale structure.
#' @details The graphon distance is computed as the L2 norm between the
#' graphons of the two structures. This function takes into account the
#' possible permutation of the blocks and returns the minimum distance.
#'
#' @keywords internal
dist_graphon_bipartite_all_permutations <- function(pis, rhos, alphas) {
  # Extracting to pi1 and pi2
  pi1 <- pis[[1]]
  pi2 <- pis[[2]]

  # Extracting to rho1 and rho2
  rho1 <- rhos[[1]]
  rho2 <- rhos[[2]]

  # Extracting to alpha1 and alpha2
  alpha1 <- alphas[[1]]
  alpha2 <- alphas[[2]]

  # Compute all possible permutations
  perms_pi1 <- gtools::permutations(length(pi1), length(pi1), seq_along(pi1))
  # Remove reversed
  # perms_pi1 <- matrix(perms_pi1[perms_pi1[, 1] < perms_pi1[, ncol(perms_pi1)], ],
  #   ncol = length(pi1)
  # )
  perms_pi2 <- gtools::permutations(length(pi2), length(pi2), seq_along(pi2))
  # Remove reversed
  # perms_pi2 <- matrix(perms_pi2[perms_pi2[, 1] < perms_pi2[, ncol(perms_pi2)], ],
  #   ncol = length(pi2)
  # )
  perms_rho1 <- gtools::permutations(length(rho1), length(rho1), seq_along(rho1))
  #  Remove reversed
  # perms_rho1 <- matrix(perms_rho1[perms_rho1[, 1] < perms_rho1[, ncol(perms_rho1)], ],
  #   ncol = length(rho1)
  # )
  perms_rho2 <- gtools::permutations(length(rho2), length(rho2), seq_along(rho2))
  # Remove reversed
  # perms_rho2 <- matrix(perms_rho2[perms_rho2[, 1] < perms_rho2[, ncol(perms_rho2)], ],
  #   ncol = length(rho2)
  # )

  dist_tensor <- array(0, dim = c(nrow(perms_pi1), nrow(perms_pi2), nrow(perms_rho1), nrow(perms_rho2)))
  dimnames(dist_tensor) <- list(
    paste0("pi1.", seq_len(nrow(perms_pi1))),
    paste0("pi2.", seq_len(nrow(perms_pi2))),
    paste0("rho1.", seq_len(nrow(perms_rho1))),
    paste0("rho2.", seq_len(nrow(perms_rho2)))
  )
  for (perm_pi1 in seq_len(nrow(perms_pi1))) {
    for (perm_pi2 in seq_len(nrow(perms_pi2))) {
      for (perm_rho1 in seq_len(nrow(perms_rho1))) {
        for (perm_rho2 in seq_len(nrow(perms_rho2))) {
          # Compute distance for current permutations
          current_distance <- graphon_distance_bipartite(
            pis = list(pi1[perms_pi1[perm_pi1, ]], pi2[perms_pi2[perm_pi2, ]]),
            rhos = list(rho1[perms_rho1[perm_rho1, ]], rho2[perms_rho2[perm_rho2, ]]),
            alphas = list(
              alpha1[perms_pi1[perm_pi1, ], perms_rho1[perm_rho1, ]],
              alpha2[perms_pi2[perm_pi2, ], perms_rho2[perm_rho2, ]]
            )
          )
          # Update minimum distance if needed
          dist_tensor[perm_pi1, perm_pi2, perm_rho1, perm_rho2] <- current_distance
        }
      }
    }
  }
  return(min(dist_tensor))
}

#' Graphon distance for bipartite SBM using identifiability
#' properties of the graphons
#' @param pis A list of two probability vectors (row)
#' @param rhos A list of two probability vectors (columns)
#' @param alphas A list of two connectivity matrices
#'
#' @return The graphon distance between two mesoscale structure.
#' @keywords internal
dist_graphon_bipartite_marginals <- function(pis, rhos, alphas) {
  # Extract pi and rho
  pi1 <- as.vector(pis[[1]])
  pi2 <- as.vector(pis[[2]])
  rho1 <- as.vector(rhos[[1]])
  rho2 <- as.vector(rhos[[2]])

  # Extract alpha
  alpha1 <- matrix(alphas[[1]], nrow = length(pi1), ncol = length(rho1))
  alpha2 <- matrix(alphas[[2]], nrow = length(pi2), ncol = length(rho2))

  # Compute order
  row_marginal1 <- as.vector(rho1 %*% t(alpha1))
  row_marginal2 <- as.vector(rho2 %*% t(alpha2))
  col_marginal1 <- as.vector(pi1 %*% alpha1)
  col_marginal2 <- as.vector(pi2 %*% alpha2)

  row_order1 <- order(row_marginal1, decreasing = TRUE)
  row_order2 <- order(row_marginal2, decreasing = TRUE)

  col_order1 <- order(col_marginal1, decreasing = TRUE)
  col_order2 <- order(col_marginal2, decreasing = TRUE)

  # # If there are ties, we need to reorder the columns too
  # if (any(rle(row_order1)[["lengths"]] > 1L)) {
  #   col_order1 <- order(col_marginal1, decreasing = TRUE)
  # } else {
  #   col_order1 <- seq_len(length(col_marginal1))
  # }
  # if (any(rle(row_order2)[["lengths"]] > 1L)) {
  #   col_order2 <- order(col_marginal2, decreasing = TRUE)
  # } else {
  #   col_order2 <- seq_len(length(col_marginal2))
  # }

  # Reorder all parameters
  pi1 <- pi1[row_order1]
  pi2 <- pi2[row_order2]
  rho1 <- rho1[col_order1]
  rho2 <- rho2[col_order2]
  alpha1 <- alpha1[row_order1, col_order1]
  alpha2 <- alpha2[row_order2, col_order2]
  return(graphon_distance_bipartite(list(pi1, pi2), list(rho1, rho2), list(alpha1, alpha2)))
}

matrix_distance_graphon_bipartite <- function(parameters_list) {
  M <- length(parameters_list)
  dist_matrix <- matrix(0, nrow = M, ncol = M)
  for (m1 in seq(1, M)) {
    for (m2 in seq(1, m1)) {
      dist_matrix[m1, m2] <- dist_graphon_bipartite_marginals(pis = list(
        parameters_list[[m1]]$pi[[1]],
        parameters_list[[m2]]$pi[[1]]
      ), rhos = list(
        parameters_list[[m1]]$rho[[1]],
        parameters_list[[m2]]$rho[[1]]
      ), alphas = list(
        parameters_list[[m1]]$alpha[[1]],
        parameters_list[[m2]]$alpha[[1]]
      ))
    }
  }
  return(t(dist_matrix) + dist_matrix)
}

#' Graphon distance for bipartite SBM using symmetrization
#'
#' @param pis A list of two probability vectors (row)
#' @param rhos A list of two probability vectors (columns)
#' @param alphas A list of two connectivity matrices
#'
#' @return The graphon distance between two mesoscale structure.
#' @keywords internal
dist_graphon_bipartite_symmetrization <- function(pis, rhos, alphas) {
  alpha1 <- alphas[[1]]
  alpha2 <- alphas[[2]]
  pi1 <- pis[[1]]
  pi2 <- pis[[2]]
  rho1 <- rhos[[1]]
  rho2 <- rhos[[2]]

  # Symmetrize alpha
  unipart_alpha1 <- matrix(0,
    nrow = length(pi1) + length(rho1),
    ncol = length(pi1) + length(rho1)
  )
  unipart_alpha1[seq_along(pi1), (length(pi1) + 1):ncol(unipart_alpha1)] <- alpha1
  unipart_alpha1[(length(pi1) + 1):nrow(unipart_alpha1), seq_along(pi1)] <- t(alpha1)
  unipart_pi1 <- c(seq_along(pi1), rep(0, length(rho1)))
  unipart_rho1 <- c(rep(0, length(pi1)), seq_along(rho1))
  unipart_prob1 <- c(pi1, rho1)
  unipart_order1 <- order(unipart_prob1 %*% unipart_alpha1, decreasing = TRUE)
  row_order1 <- unipart_pi1[unipart_order1]
  row_order1 <- row_order1[row_order1 != 0]
  col_order1 <- unipart_rho1[unipart_order1]
  col_order1 <- col_order1[col_order1 != 0]

  unipart_alpha2 <- matrix(0,
    nrow = length(pi2) + length(rho2),
    ncol = length(pi2) + length(rho2)
  )
  unipart_alpha2[seq_along(pi2), (length(pi2) + 1):ncol(unipart_alpha2)] <- alpha2
  unipart_alpha2[(length(pi2) + 1):nrow(unipart_alpha2), seq_along(pi2)] <- t(alpha2)
  unipart_pi2 <- c(seq_along(pi2), rep(0, length(rho2)))
  unipart_rho2 <- c(rep(0, length(pi2)), seq_along(rho2))
  unipart_prob2 <- c(pi2, rho2)
  unipart_order2 <- order(unipart_prob2 %*% unipart_alpha2, decreasing = TRUE)
  row_order2 <- unipart_pi2[unipart_order2]
  row_order2 <- row_order2[row_order2 != 0]
  col_order2 <- unipart_rho2[unipart_order2]
  col_order2 <- col_order2[col_order2 != 0]

  pi1 <- pi1[row_order1]
  pi2 <- pi2[row_order2]
  rho1 <- rho1[col_order1]
  rho2 <- rho2[col_order2]
  alpha1 <- alpha1[row_order1, col_order1]
  alpha2 <- alpha2[row_order2, col_order2]

  return(graphon_distance_bipartite(list(pi1, pi2), list(rho1, rho2), list(alpha1, alpha2)))
}

#' Graphon distance for unipartite SBM
#' @param pis A list of two probability vectors
#' @param alphas A list of two connectivity matrices
#'
#' @importFrom utils head tail
#'
#' @return The graphon distance between two mesoscale structure.
#' @details The graphon distance is computed as the L2 norm between the
#' graphons of the two structures. Please note that this does not take into
#' account the possible permutation of the blocks.
#' @keywords internal
graphon_distance_unipartite <- function(pis, alphas) {
  # Extracting to pi1 and pi2
  pi1 <- as.vector(pis[[1]])
  pi2 <- as.vector(pis[[2]])


  # Extracting to alpha1 and alpha2
  alpha1 <- matrix(alphas[[1]], nrow = length(pi1), ncol = length(pi1))
  alpha2 <- matrix(alphas[[2]], nrow = length(pi2), ncol = length(pi2))
  dimnames(alpha1) <- list(
    paste0("row1.", seq_len(nrow(alpha1))),
    paste0("trow1.", seq_len(ncol(alpha1)))
  )
  dimnames(alpha2) <- list(
    paste0("row2.", seq_len(nrow(alpha2))),
    paste0("trow2.", seq_len(ncol(alpha2)))
  )

  # Cumsums of pi1 pi2 rho1 rho2
  pi1_cumsum <- c(0, cumsum(pi1))
  names(pi1_cumsum) <- paste0("row1.", seq_len(length(pi1_cumsum)) - 1)
  pi2_cumsum <- c(0, cumsum(pi2))
  names(pi2_cumsum) <- paste0("row2.", seq_len(length(pi2_cumsum)) - 1)


  min_pis <- outer(tail(pi1_cumsum, -1), tail(pi2_cumsum, -1), pmin)

  max_pis <- outer(head(pi1_cumsum, -1), head(pi2_cumsum, -1), pmax)

  # All possible surfaces
  inter_surf <- outer(pmax(min_pis - max_pis, 0), t(pmax(min_pis - max_pis, 0)))
  inter_surf <- aperm(inter_surf, c(1, 4, 2, 3)) # Adjusting dimensions

  # All possible connectivity diffrences
  diff_alpha <- outer(alpha1, alpha2, "-")

  # Compute the distance
  sum(inter_surf * (diff_alpha^2))
}

#' Graphon distance for bipartite SBM over all permutations of the blocks
#' @param pis A list of two probability vectors (row)
#' @param alphas A list of two connectivity matrices
#'
#' @return The graphon distance between two mesoscale structure.
#' @details The graphon distance is computed as the L2 norm between the
#' graphons of the two structures. This function takes into account the
#' possible permutation of the blocks and returns the minimum distance.
#'
#' @keywords internal
dist_graphon_unipartite_all_permutations <- function(pis, alphas) {
  # Extracting to pi1 and pi2
  pi1 <- pis[[1]]
  pi2 <- pis[[2]]

  # Extracting to alpha1 and alpha2
  alpha1 <- alphas[[1]]
  alpha2 <- alphas[[2]]

  # Compute all possible permutations
  perms_pi1 <- gtools::permutations(length(pi1), length(pi1), seq_along(pi1))
  # Remove reversed
  # perms_pi1 <- matrix(perms_pi1[perms_pi1[, 1] < perms_pi1[, ncol(perms_pi1)], ],
  #   ncol = length(pi1)
  # )
  perms_pi2 <- gtools::permutations(length(pi2), length(pi2), seq_along(pi2))
  # Remove reversed
  # perms_pi2 <- matrix(perms_pi2[perms_pi2[, 1] < perms_pi2[, ncol(perms_pi2)], ],
  #   ncol = length(pi2)
  # )

  dist_tensor <- array(0, dim = c(nrow(perms_pi1), nrow(perms_pi2)))
  dimnames(dist_tensor) <- list(
    paste0("pi1.", seq_len(nrow(perms_pi1))),
    paste0("pi2.", seq_len(nrow(perms_pi2)))
  )
  for (perm_pi1 in seq_len(nrow(perms_pi1))) {
    for (perm_pi2 in seq_len(nrow(perms_pi2))) {
      # Compute distance for current permutations
      current_distance <- graphon_distance_unipartite(
        pis = list(pi1[perms_pi1[perm_pi1, ]], pi2[perms_pi2[perm_pi2, ]]),
        alphas = list(
          alpha1[perms_pi1[perm_pi1, ], perms_pi1[perm_pi1, ]],
          alpha2[perms_pi2[perm_pi2, ], perms_pi2[perm_pi2, ]]
        )
      )
      # Update minimum distance if needed
      dist_tensor[perm_pi1, perm_pi2] <- current_distance
    }
  }
  return(min(dist_tensor))
}

#' Graphon distance for bipartite SBM using identifiability
#' properties of the graphons
#' @param pis A list of two probability vectors
#' @param alphas A list of two connectivity matrices
#'
#' @return The graphon distance between two mesoscale structure.
#' @keywords internal
dist_graphon_unipartite_marginals <- function(pis, alphas) {
  # Extract pi and rho
  pi1 <- as.vector(pis[[1]])
  pi2 <- as.vector(pis[[2]])


  # Extract alpha
  alpha1 <- matrix(alphas[[1]], nrow = length(pi1), ncol = length(pi1))
  alpha2 <- matrix(alphas[[2]], nrow = length(pi2), ncol = length(pi2))

  # Compute order
  row_marginal1 <- as.vector(pi1 %*% t(alpha1))
  row_marginal2 <- as.vector(pi2 %*% t(alpha2))

  row_order1 <- order(row_marginal1, decreasing = TRUE)
  row_order2 <- order(row_marginal2, decreasing = TRUE)

  # Reorder all parameters
  pi1 <- pi1[row_order1]
  pi2 <- pi2[row_order2]
  alpha1 <- alpha1[row_order1, row_order1]
  alpha2 <- alpha2[row_order2, row_order2]
  return(graphon_distance_unipartite(list(pi1, pi2), list(alpha1, alpha2)))
}

#' Generate a unipartite network
#'
#' @param n the number of nodes
#' @param pi a vector of probability to belong to the clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param distribution the emission distribution, either "bernoulli" or
#' "poisson"
#' @param return_memberships Boolean, should return memberships or not.
#'
#' @return An adjacency matrix

#' @noMd
#' @noRd
generate_unipartite_network <- function(
    n, pi, alpha,
    distribution = "bernoulli",
    return_memberships = FALSE) {
  stopifnot(
    "All alpha coefficients must be positive" = all(alpha >= 0),
    "With bernoulli, the alpha must be between 0 and 1" = (
      distribution == "poisson" |
        (distribution == "bernoulli" & all(alpha >= 0L & alpha <= 1L))),
    "All pi must be between 0 and 1" = all(pi >= 0L & pi <= 1L),
    "Pi must sum to one" = all.equal(sum(pi), 1L)
  )
  cluster_memberships <- rmultinom(n, size = 1, prob = pi)
  node_node_interaction_parameter <- t(cluster_memberships) %*% alpha %*% cluster_memberships

  #  Here we switch on the distributions
  adjacency_matrix <- matrix(
    switch(distribution,
      "bernoulli" = {
        rbinom(length(node_node_interaction_parameter),
          size = 1, prob = node_node_interaction_parameter
        )
      },
      "poisson" = {
        rpois(length(node_node_interaction_parameter),
          lambda = node_node_interaction_parameter
        )
      },
      stop("distribution must be one of either 'bernoulli' or 'poisson'")
    ),
    nrow = nrow(node_node_interaction_parameter)
  )

  if (return_memberships) {
    return(list(
      adjacency_matrix = adjacency_matrix,
      block_memberships = cluster_memberships
    ))
  } else {
    return(adjacency_matrix)
  }
}

#' Generate collection of unipartite
#'
#' @param n the number of nodes or a vector of the nodes per network
#' @param pi a vector of probability to belong to the clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param M the number of networks to generate
#' @param distribution the emission distribution, either "bernoulli" or
#' "poisson"
#' @param return_memberships Boolean, should return memberships or not.
#'
#' @details If n is a single value, this value will be replicated for each of
#' the M networks. If it is a vector, it must be of size M, specifying the
#' number of nodes for each network.
#'
#' @return A list of M lists, which contains : $adjacency_matrix, $clustering
#'
#' @export
generate_unipartite_collection <- function(
    n, pi, alpha, M,
    distribution = "bernoulli",
    return_memberships = FALSE) {
  if (length(n) == 1) {
    n <- rep(n, M)
  }

  # Check if n is the correct length
  if (length(n) != M) {
    stop(
      "The length of n is not correct ! It should be : ",
      M, " values and it is ", length(n)
    )
  }
  # Generate the networks
  out <- lapply(seq.int(M), function(m) {
    generate_unipartite_network(
      n = n[[m]],
      pi = pi,
      alpha = alpha,
      distribution = distribution
    )
  })
  return(out)
}

#' Generate a bipartite network
#'
#' @param nr the number of row nodes
#' @param nc the number of col nodes
#' @param pi a vector of probability to belong to the row clusters
#' @param rho a vector of probability to belong to the columns clusters
#' @param alpha the matrix of connectivity between two clusters
#'
#' @return An incidence matrix
#'
#' @noMd
#' @noRd
generate_bipartite_network <- function(
    nr, nc, pi, rho, alpha, distribution = "bernoulli",
    return_memberships = FALSE) {
  stopifnot(
    "All alpha coefficients must be positive" = all(alpha >= 0L),
    "With bernoulli, the alpha must be between 0 and 1" = (
      distribution == "poisson" |
        (distribution == "bernoulli" & all(alpha >= 0L & alpha <= 1L))),
    "All pi must be between 0 and 1" = all(pi >= 0L & pi <= 1L),
    "Pi must sum to one" = all.equal(sum(pi), 1L),
    "All rho must be between 0 and 1" = all(rho >= 0L & rho <= 1L),
    "Rho must sum to one" = all.equal(sum(rho), 1L)
  )

  rowblocks_memberships <- rmultinom(nr, size = 1, prob = pi)
  colblocks_memberships <- rmultinom(nc, size = 1, prob = rho)
  node_node_interaction_parameter <- t(rowblocks_memberships) %*%
    alpha %*%
    colblocks_memberships
  incidence_matrix <- matrix(
    switch(distribution,
      "poisson" = {
        rpois(length(node_node_interaction_parameter),
          lambda = node_node_interaction_parameter
        )
      },
      "bernoulli" = {
        rbinom(length(node_node_interaction_parameter),
          size = 1, prob = node_node_interaction_parameter
        )
      },
      stop("distribution must be one of either 'bernoulli' or 'poisson'")
    ),
    nrow = nrow(node_node_interaction_parameter)
  )

  if (return_memberships) {
    return(list(
      incidence_matrix = incidence_matrix,
      row_blockmemberships = as.vector(
        c(seq.int(length(pi))) %*% rowblocks_memberships
      ), # We reverse the one hot encoding
      col_blockmemberships = as.vector(
        c(seq.int(length(rho))) %*% colblocks_memberships
      ) # We reverse the one hot encoding
    ))
  } else {
    return(incidence_matrix)
  }
}

#' Generate collection of bipartite networks
#'
#' @param nr the number of row nodes  or a vector specifying the
#' number of row nodes for each of the M networks
#' @param nc the number of column nodes  or a vector specifying the
#' number of column nodes for each of the M networks
#' @param pi a vector of probability to belong to the row clusters
#' @param rho a vector of probability to belong to the columns clusters
#' @param alpha the matrix of connectivity between two clusters
#' @param M the number of networks to generate
#' @param model the colBiSBM model to use. Available: "iid", "pi", "rho",
#' "pirho"
#' @param distribution the emission distribution : "bernoulli" or "poisson"
#' @param return_memberships a boolean which choose whether the function returns
#' a list containing the memberships and the incidence matrices or just the
#' incidence matrices. Defaults to FALSE, only the matrices are returned.
#'
#' @details the model parameters if set to any other than iid will shuffle the
#' provided pi and rho
#'
#' @return A list of M lists, which contains : $incidence_matrix, $row_blockmemberships, $col_blockmemberships
#'
#' @export
generate_bipartite_collection <- function(
    nr, nc, pi, rho, alpha, M,
    model = "iid",
    distribution = "bernoulli",
    return_memberships = FALSE) {
  out <- list()

  # Check if nr and nc are vectors
  if (length(nr) == 1) {
    nr <- rep(nr, M)
  }
  if (length(nc) == 1) {
    nc <- rep(nc, M)
  }

  # Check if nr and nc are the correct length
  if (length(nr) != M) {
    stop(
      "The length of nr is not correct ! It should be : ",
      M, " values and it is ", length(nr)
    )
  }
  if (length(nc) != M) {
    stop(
      "The length of nc is not correct ! It should be : ",
      M, " values and it is ", length(nc)
    )
  }

  switch(model,
    "iid" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = pi,
          rho = rho,
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "pi" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = sample(pi),
          rho = rho,
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "rho" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = pi,
          rho = sample(rho),
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    "pirho" = {
      out <- lapply(seq.int(M), function(m) {
        generate_bipartite_network(
          nr = nr[[m]],
          nc = nc[[m]],
          pi = sample(pi),
          rho = sample(rho),
          alpha = alpha,
          distribution = distribution,
          return_memberships = return_memberships
        )
      })
    },
    stop("Error unknown model. Must be one of : iid, pi, rho, pirho.")
  )

  return(out)
}

#' Perform a spectral clustering
#'
#' @importFrom stats kmeans
#'
#' @param X an Adjacency matrix
#' @param K the number of clusters
#'
#' @noMd
#' @noRd
#'
#' @return A vector : The clusters labels
spectral_clustering <- function(X, K) {
  if (K == 1) {
    return(rep(1L, nrow(X)))
  }
  n <- nrow(X)
  if (n < 3) {
    return(rep(1, n))
  }
  X[X == -1] <- NA
  isolated <- which(rowSums(X, na.rm = TRUE) == 0)
  connected <- setdiff(seq(n), isolated)
  X <- X[connected, connected]
  if (!isSymmetric(X)) {
    X <- 1 * ((X + t(X)) > 0) # .5 * (X + t(X))
  }
  D_moins1_2 <- diag(1 / sqrt(rowSums(X, na.rm = TRUE) + 1))
  X[is.na(X)] <- mean(X, na.rm = TRUE)
  Labs <- D_moins1_2 %*% X %*% D_moins1_2
  specabs <- eigen(Labs, symmetric = TRUE)
  if (K >= nrow(X)) {
    message("Too many clusters for Spectral Clustering")
    K_old <- K
    K <- nrow(X) - 1
  }
  if (K >= 2) {
    index <- rev(order(abs(specabs$values)))[1:K]
    U <- specabs$vectors[, index]
    U <- U / rowSums(U**2)**(1 / 2)
    U[is.na(U)] <- 0
    cl <- stats::kmeans(U, K, iter.max = 100, nstart = 100)$cluster
  } else {
    cl <- rep(1, nrow(X))
  }
  index <- rev(order(abs(specabs$values)))[1:K]
  U <- specabs$vectors[, index]
  U <- U / rowSums(U**2)**(1 / 2)
  U[is.na(U)] <- 0
  U[is.nan(U)] <- 0
  U[is.infinite(U)] <- 0
  print(U)
  cl <- stats::kmeans(U[!is.na(U)], K, iter.max = 100, nstart = 100)$cluster
  clustering <- rep(1, n)
  clustering[connected] <- cl
  clustering[isolated] <- which.min(rowsum(rowSums(X, na.rm = TRUE), cl))
  clustering[isolated] <- which.min(rowsum(rowSums(X, na.rm = TRUE), cl))
  return(clustering)
}



#' Perform a spectral bi-clustering, clusters by row
#' and by columns independently
#'
#' Relies on the spectral_clustering function defined above
#'
#' @param X an Incidence matrix
#' @param K the two numbers of clusters
#'
#' @noMd
#' @noRd
#'
#' @return A list of two vectors : The clusters labels.
#' They are accessed using $row_clustering and $col_clustering
spectral_biclustering <- function(X, K) {
  # Trivial clustering : everyone is part of the cluster
  if (all(K == c(1, 1))) {
    return(list(
      row_clustering = rep(1, nrow(X)),
      col_clustering = rep(1, ncol(X))
    ))
  }

  # Extracts the number of clusters
  K1 <- K[1] # Row clusters
  K2 <- K[2] # Column clusters

  row_adjacency_matrix <- tcrossprod(X)
  row_clustering <- spectral_clustering(row_adjacency_matrix, K1)

  col_adjacency_matrix <- crossprod(X)
  col_clustering <- spectral_clustering(col_adjacency_matrix, K2)

  return(list(row_clustering = row_clustering, col_clustering = col_clustering))
}

# TODO : Modify the algorithm to use the rectangular matrix and its transpose
bipartite_hierarchic_clustering <- function(X, K) {
  # Trivial clustering : everyone is part of the cluster
  if (all(K == c(1, 1))) {
    return(list(
      row_clustering = rep(1, nrow(X)),
      col_clustering = rep(1, ncol(X))
    ))
  }

  # Extracts the number of clusters
  K1 <- K[1] # Row clusters
  K2 <- K[2] # Column clusters

  row_adjacency_matrix <- tcrossprod(X)
  row_clustering <- hierarClust(X, K1)

  col_adjacency_matrix <- crossprod(X)
  col_clustering <- hierarClust(t(X), K2)

  return(list(row_clustering = row_clustering, col_clustering = col_clustering))
}

#' Perform a Hierarchical Clustering
#' @importFrom stats cutree dist hclust
#' @importFrom ape additive
#' @param X An Adjacency Matrix
#' @param K the number of wanted clusters
#'
#' @noMd
#' @noRd
#'
#' @return A vector : The clusters labels
hierarClust <- function(X, K) {
  if (K == 1) {
    return(rep(1L, nrow(X)))
  }
  # distance <- stats::dist(x = X, method = "manhattan")
  # X[X == -1] <- NA
  # distance[which(A == 1)] <- distance[which(A == 1)] - 2
  # distance <- stats::as.dist(ape::additive(distance))
  diss <- cluster::daisy(x = X, metric = "manhattan", warnBin = FALSE)
  if (!any(is.na(diss))) {
    clust <- cluster::agnes(x = X, metric = "manhattan", method = "ward")
  } else {
    return(rep(1L, nrow(X)))
  }
  # clust    <- stats::hclust(d = distance , method = "ward.D2")
  return(stats::cutree(tree = clust, k = K))
}

#' Split a list of clusters
#'
#' @param X an adjacency matrix
#' @param Z a vector of cluster memberships
#' @param Q The number of maximal clusters
#'
#' @noMd
#' @noRd
#'
#' @return A list of Q clustering of Q+1 clusters
split_clust <- function(X, Z, Q, is_bipartite = FALSE) {
  Z_split <- lapply(
    seq(Q),
    FUN = function(q) {
      if (sum(Z == q) <= 3) {
        return(Z)
      }
      Z_new <- Z
      mx <- ifelse(is_bipartite, mean(X[Z == q, ], na.rm = TRUE), mean(X[Z == q, Z == q], na.rm = TRUE))
      X[is.na(X)] <- mx
      if (isSymmetric.matrix(X) || is_bipartite) {
        Xsub <- X[Z == q, , drop = FALSE]
      } else {
        Xsub <- cbind(X[Z == q, , drop = FALSE], t(X[, Z == q, drop = FALSE]))
      }
      if (nrow(unique(Xsub, MARGIN = 1)) <= 3) {
        return(Z)
      }
      C <- stats::kmeans(x = Xsub, centers = 2)$cluster
      # C  <-  hierarClust(Xsub, 2)# + Q#stats::kmeans(x = .5 * (X[Z==q,] + t(X[,Z==q])), centers = 2)$cluster + Q
      if (length(unique(C)) == 2) {
        if (!is_bipartite) {
          p1 <- c(
            mean(X[Z == q, , drop = FALSE][C == 1, ]),
            mean(X[, Z == q, drop = FALSE][, C == 1])
          )
          p2 <- c(
            mean(X[Z == q, , drop = FALSE][C == 2, ]),
            mean(X[, Z == q, drop = FALSE][, C == 2])
          )
        } else {
          p1 <- c(
            mean(X[Z == q, , drop = FALSE][C == 1, ])
          )
          p2 <- c(
            mean(X[Z == q, , drop = FALSE][C == 2, ])
          )
        }
        c <- which.max(abs(p1 - p2))
        md <- sample(
          x = 2, size = 2, replace = FALSE,
          prob = c(
            max(1 / ncol(Xsub), p1[c]),
            max(1 / ncol(Xsub), p2[c])
          )
        )
        Z_new[Z == q][C == which.min(md)] <- Q + 1
      }
      # Z_new[Z==q]  <-  stats::kmeans(x = Xsub, centers = 2)$cluster + Q
      # Z_new[Z_new == Q + 2]  <-  q
      return(Z_new)
    }
  )
  Z_split <- Z_split[which(sapply(X = Z_split, FUN = function(x) !is.null(x)))]
  return(Z_split)
}

#' Merge a list of clusters
#'
#' @importFrom utils combn
#'
#' @param Z a vector of cluster memberships
#' @param Q the number of original clusters
#'
#' @return A list of Q(Q-1)/2 clustering of Q-1 clusters
#'
#'  @noMd
#' @noRd
merge_clust <- function(Z, Q) {
  Z_merge <- lapply(
    X = 1:choose(Q, 2),
    FUN = function(q) {
      Z[Z == utils::combn(Q, 2)[2, q]] <- utils::combn(Q, 2)[1, q]
      if (utils::combn(Q, 2)[2, q] < Q) {
        Z[Z > utils::combn(Q, 2)[2, q]] <- Z[Z > utils::combn(Q, 2)[2, q]] - 1
      }
      return(Z)
    }
  )
  return(Z_merge)
}



#
# Fonction interne a VEM
#

F.bern <- function(X, alpha, tau) {
  return(X %*% tau %*% log(alpha) +
    (1 - X - diag(1, nrow(X))) %*% tau %*% log(1 - alpha))
}

rotate <- function(x) t(apply(x, 2, rev))

#
dist_param <- function(param, param_old) {
  sqrt(sum((param - param_old)**2))
}

#' Title
#'
#' @param X An adjacency matrix
#' @param K An integer, the number of folds
#'
#' @return A matrix of the same size than X with class integer as coefficient
#'
#' @noMd
#' @noRd
build_fold_matrix <- function(X, K) {
  n <- ncol(X)
  arrange <- sample.int(n)
  labels <- cut(seq(n), breaks = K, labels = FALSE)
  fold_matrix <- diag(n)
  for (i in seq(n)) {
    fold_matrix[i, ] <- (labels + labels[i]) %% K
  }
  fold_matrix <- fold_matrix[arrange, arrange] + 1
  diag(fold_matrix) <- 0
  return(fold_matrix)
}


.xlogx <- function(x) {
  ifelse(x < 2 * .Machine$double.eps, 0, x * log(x))
}
.xlogy <- function(x, y, eps = NULL) {
  ifelse(x < 2 * .Machine$double.eps, 0, x * .log(y, eps = eps))
}
.quadform <- function(x, y) tcrossprod(x %*% y, x)
.tquadform <- function(x, y) crossprod(x, y %*% x)
logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))
.logit <- function(x, eps = NULL) {
  if (is.null(eps)) {
    res <- log(x / (1 - x))
  } else {
    res <- log(pmax(pmin(x, 1 - eps), eps) / pmax(pmin(1 - x, 1 - eps), eps))
  }
  return(res)
}


.threshold <- function(x, eps = 1e-9) {
  #  x <- .softmax(x)
  x[x < eps] <- eps
  x[x > 1 - eps] <- 1 - eps
  x <- x / .rowSums(x, nrow(x), ncol(x))
  x
}

.softmax <- function(x) {
  # x_max <- apply(x, 1, max)
  # OPTIM
  x_max <- matrixStats::rowMaxs(x)
  x <- exp(x - x_max)
  x <- x / .rowSums(x, nrow(x), ncol(x))
  x
}
.log <- function(x, eps = NULL) {
  if (is.null(eps)) {
    res <- log(x)
  } else {
    # OPTIM
    # x[x >= 1 - eps] <- 1 - eps
    x[x <= eps] <- eps
    res <- log(x)
  }
  return(res)
}
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x), Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}

.rev_one_hot <- function(X) {
  return(as.vector(max.col(X)))
}

#' Reorder colBiSBM parameters
#'
#' @param model a fitBipartiteSBMPop object on which the reordering is performed
#'
#' @export
reorder_parameters <- function(model) {
  Z_label_switch <- function(Z, new_order) {
    # Create a mapping of old labels to new labels
    old_names <- names(Z)
    label_map <- setNames(new_order, unique(Z))

    # Use the mapping to replace labels in the vector
    switched_labels <- label_map[Z]
    names(switched_labels) <- old_names
    return(switched_labels)
  }
  out_model <- model$clone()
  if (all(out_model$Q == c(1, 1))) {
    return(out_model)
  }
  mean_pi <- sapply(out_model$pi, function(pi) pi[[1]])
  if (out_model$Q[1] > 1) {
    mean_pi <- matrixStats::rowMeans2(mean_pi)
  } else {
    mean_pi <- 1
  }
  mean_rho <- sapply(out_model$pi, function(pi) pi[[2]])
  if (out_model$Q[2] > 1) {
    mean_rho <- matrixStats::rowMeans2(mean_rho)
  } else {
    mean_rho <- 1
  }
  # The row clustering are reordered according to their marginal distribution
  prob1 <- as.vector(mean_rho %*% t(out_model$MAP$alpha))
  p1 <- order(prob1, decreasing = TRUE)

  # The col clustering are reordered according to their marginal distribution
  prob2 <- as.vector(mean_pi %*% out_model$MAP$alpha)
  p2 <- order(prob2, decreasing = TRUE)

  # m independent
  out_model$MAP$alpha <- out_model$MAP$alpha[p1, p2, drop = FALSE]
  out_model$Calpha <- out_model$Calpha[p1, p2, drop = FALSE]
  out_model$alpha <- out_model$alpha[p1, p2, drop = FALSE]

  # m dependent
  lapply(seq.int(out_model$M), function(m) {
    # Reordering the parameters
    out_model$Cpi[[1]][, m] <- out_model$Cpi[[1]][p1, m, drop = FALSE]
    out_model$Cpi[[2]][, m] <- out_model$Cpi[[2]][p2, m, drop = FALSE]

    out_model$pim[[m]][[1]] <- out_model$pim[[m]][[1]][p1, drop = FALSE]
    out_model$pim[[m]][[2]] <- out_model$pim[[m]][[2]][p2, drop = FALSE]

    out_model$pi[[m]][[1]] <- out_model$pi[[m]][[1]][p1, drop = FALSE]
    out_model$pi[[m]][[2]] <- out_model$pi[[m]][[2]][p2, drop = FALSE]

    out_model$emqr[m, , ] <- out_model$emqr[m, p1, p2, drop = FALSE]
    out_model$nmqr[m, , ] <- out_model$nmqr[m, p1, p2, drop = FALSE]
    out_model$alpham[[m]] <- matrix(out_model$alpham[[m]], out_model$Q[1], out_model$Q[2])[p1, p2, drop = FALSE]
    out_model$tau[[m]][[1]] <- out_model$tau[[m]][[1]][, p1, drop = FALSE]
    out_model$tau[[m]][[2]] <- out_model$tau[[m]][[2]][, p2, drop = FALSE]
    out_model$Z[[m]][[1]] <- Z_label_switch(out_model$Z[[m]][[1]], p1)
    out_model$Z[[m]][[2]] <- Z_label_switch(out_model$Z[[m]][[2]], p2)

    # MAP parameters
    # Work needed to relabel correctly!
    out_model$MAP$Z[[m]][[1]] <- Z_label_switch(out_model$MAP$Z[[m]][[1]], p1)
    out_model$MAP$Z[[m]][[2]] <- Z_label_switch(out_model$MAP$Z[[m]][[2]], p2)
    out_model$MAP$emqr[m, , ] <- out_model$MAP$emqr[m, p1, p2, drop = FALSE]
    out_model$MAP$nmqr[m, , ] <- out_model$MAP$nmqr[m, p1, p2, drop = FALSE]
    out_model$MAP$alpham[[m]] <- matrix(out_model$MAP$alpham[[m]], out_model$Q[1], out_model$Q[2])[p1, p2, drop = FALSE]
    out_model$MAP$pim[[m]][[1]] <- out_model$MAP$pim[[m]][[1]][p1, drop = FALSE]
    out_model$MAP$pim[[m]][[2]] <- out_model$MAP$pim[[m]][[2]][p2, drop = FALSE]
    out_model$MAP$pi[[m]][[1]] <- out_model$MAP$pi[[m]][[1]][p1, drop = FALSE]
    out_model$MAP$pi[[m]][[2]] <- out_model$MAP$pi[[m]][[2]][p2, drop = FALSE]
  })
  return(out_model)
}

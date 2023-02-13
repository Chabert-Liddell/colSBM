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
  index <- rev(order(abs(specabs$values)))[1:K]
  U <- specabs$vectors[, index]
  U <- U / rowSums(U**2)**(1 / 2)
  U[is.na(U)] <- 0
  cl <- stats::kmeans(U, K, iter.max = 100, nstart = 100)$cluster
  clustering <- rep(1, n)
  clustering[connected] <- cl
  clustering[isolated] <- which.min(rowsum(rowSums(X, na.rm = TRUE), cl))
  return(clustering)
}



# FIXME : implement a spectral bi-clustering, check if it works consistently with co-clustering
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
#' @return A list of two vectors : The clusters labels
spectral_biclustering <- function(X, K) {
  if (K == c(1, 1)) {
    return(list(row_clustering = rep(1, nrow(X)), col_clustering = rep(1, ncol(X))))
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

# TODO : implement a spectral co-clustering
#' Perform a spectral co-clustering, clusters by row
#' and by columns
#'
#' @importFrom stats kmeans
#'
#' @param X an Incidence matrix
#' @param K the two numbers of clusters
#'
#' @noMd
#' @noRd
#'
#' @return A list of two vectors : The clusters labels
spectral_coclustering <- function(X, K) {
  # placeholder
}

# TODO : implement CAH bi-clustering (GREMLINS)
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

#' Merge a list of clusters
#'
#' @param X an adjacency matrix
#' @param Z a vector of cluster memberships
#' @param Q The number of maximal clusters
#'
#' @noMd
#' @noRd
#'
#' @return A list of Q clustering of Q+1 clusters
split_clust <- function(X, Z, Q) {
  Z_split <- lapply(
    X = seq(Q),
    FUN = function(q) {
      if (sum(Z == q) <= 3) {
        return(Z)
      }
      Z_new <- Z
      mx <- mean(X[Z == q, Z == q], na.rm = TRUE)
      X[is.na(X)] <- mx
      if (isSymmetric.matrix(X)) {
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
        p1 <- c(
          mean(X[Z == q, , drop = FALSE][C == 1, ]),
          mean(X[, Z == q, drop = FALSE][, C == 1])
        )
        p2 <- c(
          mean(X[Z == q, , drop = FALSE][C == 2, ]),
          mean(X[, Z == q, drop = FALSE][, C == 2])
        )
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
generate_unipartite_collection <- function(nodesPerClass, numberOfClass, numberOfNetworks, directed = F) {
  A <- list()
  for (m in seq(numberOfNetworks)) {
    n <- nodesPerClass * numberOfClass

    Z <- diag(numberOfClass) %x% matrix(1, nodesPerClass, 1)
    P <- matrix(runif(numberOfClass * numberOfClass), numberOfClass, numberOfClass)
    if (directed) {
      P[lower.tri(P)] <- t(P)[lower.tri(P)]
    }
    A[[m]] <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z)) ## adjacency matrix
    if (directed) {
      A[lower.tri(A)] <- t(A)[lower.tri(A)]
    }
  }
  A
}

generate_bipartite_collection <- function(coupleOfNodesPerClass, coupleOfNumberOfClasses, numberOfNetworks) {
  out <- list()
  out$n <- list()
  out$A <- list()
  out$Z1 <- list()
  out$Z2 <- list()
  out$P <- list()
  for (m in seq(numberOfNetworks)) {
    n <- coupleOfNodesPerClass * coupleOfNumberOfClasses # nodes
    Z1 <- diag(coupleOfNumberOfClasses[1]) %x% matrix(1, coupleOfNodesPerClass[1], 1)
    Z2 <- diag(coupleOfNumberOfClasses[2]) %x% matrix(1, coupleOfNodesPerClass[2], 1)
    P <- matrix(runif(coupleOfNumberOfClasses[1] * coupleOfNumberOfClasses[2]), coupleOfNumberOfClasses[1], coupleOfNumberOfClasses[2])
    out$A[[m]] <- 1 * (matrix(runif(n[1] * n[2]), n[1], n[2]) < Z1 %*% P %*% t(Z2)) ## adjacency matrix

    out$n[[m]] <- n
    out$Z1[[m]] <- Z1
    out$Z2[[m]] <- Z2
    out$P[[m]] <- P
  }
  out
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


.xlogx <- function(x) ifelse(x < 2 * .Machine$double.eps, 0, x * log(x))
.xlogy <- function(x, y, eps = NULL) {
  ifelse(x < 2 * .Machine$double.eps, 0, x * .log(y, eps = eps))
}
.quadform <- function(x, y) tcrossprod(x %*% y, x)
.tquadform <- function(x, y) crossprod(x, y %*% x) # TODO : check if naming is consistent with R functions
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
  x_max <- apply(x, 1, max)
  x <- exp(x - x_max)
  x <- x / .rowSums(x, nrow(x), ncol(x))
  x
}
.log <- function(x, eps = NULL) {
  if (is.null(eps)) {
    res <- log(x)
  } else {
    res <- log(pmax(pmin(x, 1 - eps), eps))
  }
  return(res)
}
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x), Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}

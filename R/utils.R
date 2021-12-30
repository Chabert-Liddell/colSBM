#' Perform a spectral clustering
#'
#' @importFrom stats kmeans
#'
#' @param X an Adjacency matrix
#' @param K the number of clusters
#'
#' @return A vector : The clusters labels
spectral_clustering <- function(X, K){
  if (K == 1) return (rep(1L, nrow(X)))
  n <- nrow(X)
  X[X == -1] <- NA
  isolated <- which(rowSums(X, na.rm = TRUE) == 0)
  connected <- setdiff(seq(n), isolated)
  X <- X[connected, connected]
  if (! isSymmetric(X)) {
    X <- 1*((X+t(X)) >0)#.5 * (X + t(X))
  }
  D_moins1_2 <- diag(1/sqrt(rowSums(X, na.rm = TRUE) + 1))
  X[is.na(X)] <- mean(X, na.rm=TRUE)
  Labs <- D_moins1_2 %*% X %*% D_moins1_2
  specabs <- eigen(Labs, symmetric = TRUE)
  index <- rev(order(abs(specabs$values)))[1:K]
  U <- specabs$vectors[,index]
  U <- U / rowSums(U**2)**(1/2)
  U[is.na(U)] <- 0
  cl <- stats::kmeans(U, K, iter.max = 100, nstart=100)$cluster
  clustering <- rep(1, n)
  clustering[connected] <- cl
  clustering[isolated] <-   which.min(rowsum(rowSums(X, na.rm = TRUE),cl))
  return(clustering)
}

#' Perform a Hierarchical Clustering
#' @importFrom stats cutree dist hclust
#' @importFrom ape additive
#' @param X An Adjacency Matrix
#' @param K the number of wanted clusters
#'
#' @return A vector : The clusters labels
hierarClust <- function(X, K){
  if (K == 1) return (rep(1L, nrow(X)))
  # distance <- stats::dist(x = X, method = "manhattan")
  # X[X == -1] <- NA
  # distance[which(A == 1)] <- distance[which(A == 1)] - 2
  # distance <- stats::as.dist(ape::additive(distance))
  diss <- cluster::daisy(x = X, metric = "manhattan")
  if (! any(is.na(diss))) {
    clust <- cluster::agnes(x = X, metric = "manhattan", method = "ward")
  } else {
    return (rep(1L, nrow(X)))
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
#' @return A list of Q clustering of Q+1 clusters
split_clust <- function(X, Z, Q) {
  Z_split <-  lapply(
    X = seq(Q),
    FUN =  function(q) {
      if (sum(Z==q) <= 3) return(Z)
      Z_new        <-  Z
 #     Z_new[Z==q]  <-  hierarClust(X[Z==q,] + t(X[,Z==q]), 2) + Q#stats::kmeans(x = .5 * (X[Z==q,] + t(X[,Z==q])), centers = 2)$cluster + Q
      mx <- mean(X[Z==q, Z==q], na.rm = TRUE)
      X[is.na(X)] <- mx
      Xsub <- cbind(X[Z==q,], t(X[,Z==q]))
      if(nrow(unique(Xsub, MARGIN =1)) <=3 ) return(Z)
      C <- stats::kmeans(x = Xsub, centers = 2)$cluster
      if (length(unique(C)) == 2) {
        md <- sample(x = 2,size =  2, replace = FALSE,
                     prob = c(max(1/ncol(Xsub), mean(Xsub[C == 1,])),
                              max(1/ncol(Xsub),mean(Xsub[C == 2,]))))
        Z_new[Z == q][C == which.min(md)] <- Q+1
      }
      # Z_new[Z==q]  <-  stats::kmeans(x = Xsub, centers = 2)$cluster + Q
      # Z_new[Z_new == Q + 2]  <-  q
      return(Z_new)
    })
  Z_split  <-  Z_split[which(sapply(X = Z_split, FUN = function(x) ! is.null(x)))]
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
merge_clust <- function(Z, Q) {
  Z_merge <- lapply(
    X = 1:choose(Q,2),
    FUN = function(q) {
      Z[Z == utils::combn(Q, 2)[2, q]] <- utils::combn(Q, 2)[1, q]
      if(utils::combn(Q, 2)[2, q] < Q){
        Z[Z > utils::combn(Q, 2)[2, q]] <- Z[Z > utils::combn(Q, 2)[2, q]] - 1
      }
      return(Z)
    }
  )
  return(Z_merge)
}



#
#Fonction interne a VEM
#

F.bern <- function(X, alpha, tau){
  return(X %*% tau %*% log(alpha) +
           (1 - X - diag(1, nrow(X))) %*% tau %*% log(1-alpha))
}

rotate <- function(x) t(apply(x, 2, rev))

#
dist_param <- function(param, param_old) {
  sqrt(sum((param-param_old)**2))
}

#' Title
#'
#' @param X An adjacency matrix
#' @param K An integer, the number of folds
#'
#' @return A matrix of the same size than X with class integer as coefficient
build_fold_matrix <- function(X, K) {
  n <- ncol(X)
  arrange     <- sample(x = seq(n))
  labels      <- cut(seq(n), breaks = K, labels = FALSE)
  fold_matrix <- diag(n)
  for (i in seq(n)) {
    fold_matrix[i, ] <- (labels + labels[i]) %% K
  }
  fold_matrix <- fold_matrix[arrange, arrange] + 1
  diag(fold_matrix) <- 0
  return(fold_matrix)
}


.xlogx      <- function(x) ifelse(x < 2*.Machine$double.eps, 0, x*log(x))
.xlogy     <- function(x, y, eps = NULL) {
  ifelse(x < 2*.Machine$double.eps, 0, x*.log(y, eps = eps))
}
.quadform  <- function(x, y)  tcrossprod(x %*% y, x)
.tquadform  <- function(x, y)  crossprod(x, y %*% x)
logistic   <- function(x) 1/(1 + exp(-x))
logit      <- function(x) log(x/(1 - x))
.logit <- function(x, eps = NULL) {
  if(is.null(eps)) {
    res <- log(x/(1-x))
  } else {
    res <- log(pmax(pmin(x, 1-eps), eps)/pmax(pmin(1-x, 1-eps), eps))
  }
  return (res)
}


.threshold <- function(x, eps = 1e-9) {
#  x <- .softmax(x)
  x[x < eps] <- eps
  x[x > 1-eps] <- 1-eps
  x <- x/.rowSums(x, nrow(x), ncol(x))
  x
}

.softmax <- function(x) {
  x_max <- apply(x, 1, max)
  x <- exp(x - x_max)
  x <- x/.rowSums(x, nrow(x), ncol(x))
  x
}
.log <- function(x, eps = NULL) {
  if(is.null(eps)) {
    res <- log(x)
  } else {
    res <- log(pmax(pmin(x, 1-eps), eps))
  }
  return (res)
}
.one_hot <- function(x, Q) {
  O <- matrix(0, length(x),Q)
  O[cbind(seq.int(length(x)), x)] <- 1
  return(O)
}

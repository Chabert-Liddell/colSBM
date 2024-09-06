# Testing the generate_bipartite_collection
set.seed(1234)

eps <- 0.05

pir <- c(0.1, 0.5, 0.2, 0.2)
pic <- c(1 / 3, 1 / 3, 1 / 3)

alpha <- matrix(c(
  0.7, 0.6, 0.6,
  0.6, 0.5, eps,
  0.6, eps, eps,
  eps, eps, 0.9
), nrow = 4, byrow = TRUE)

M <- 5

# With one nr and one nc
nr <- 120
nc <- 240

firstCol <- colSBM::generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("nc and nr are extended to match M", {
  expect_equal(length(firstCol), M)
  expect_equal(unique(lapply(firstCol, nrow)), list(nr))
  expect_equal(unique(lapply(firstCol, ncol)), list(nc))
})

# With increasing nr and nc
nr <- c(5, 10, 20, 40, 80)
nc <- 2 * nr

increasingCol <- colSBM::generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("Generate bipartite collection with vectors for nc and nr", {
  expect_equal(length(increasingCol), M)
  expect_equal(lapply(increasingCol, nrow), as.list(nr))
  expect_equal(lapply(increasingCol, ncol), as.list(nc))
  #  Wrong M
  expect_error(colSBM::generate_bipartite_collection(
    nr, nc, pir, pic, alpha,
    M + 1
  ))
  expect_error(colSBM::generate_bipartite_collection(
    c(nr, 1), nc, pir, pic, alpha,
    M
  ))
  expect_error(colSBM::generate_bipartite_collection(
    nr, c(nc, 1), pir, pic, alpha,
    M
  ))
})

test_that("Testing the different models for generating bipartite networks", {
  expect_error(
    colSBM::generate_bipartite_collection(
      nr, nc, pir, pic, alpha,
      M,
      model = "iidi"
    )
  )
  expect_no_error(
    colSBM::generate_bipartite_collection(
      nr, nc, pir, pic, alpha,
      M,
      model = "iid"
    )
  )
  expect_no_error(
    colSBM::generate_bipartite_collection(
      nr, nc, pir, pic, alpha,
      M,
      model = "pi"
    )
  )
  expect_no_error(
    colSBM::generate_bipartite_collection(
      nr, nc, pir, pic, alpha,
      M,
      model = "rho"
    )
  )
  expect_no_error(
    colSBM::generate_bipartite_collection(
      nr, nc, pir, pic, alpha,
      M,
      model = "pirho"
    )
  )
})

test_that("Testing various wrong arguments for generating bipartite networks", {
  nr <- 10L
  nc <- 10L

  alpha <- matrix(c(1L, 1L, 1L, 1L), nrow = 2)
  alpha_3 <- matrix(c(
    0.9, 0.4, 0.1,
    0.4, 0.1, 0.05,
    0.1, 0.05, 0.05
  ), byrow = TRUE, ncol = 3)
  wrong_alpha_bernoulli <- matrix(2L)
  wrong_alpha_negative <- matrix(-5L)

  pi <- c(0.5, 0.5)
  pi_sum_float <- c(1 / 3 + 0.28, 1 / 3, 1 / 3 - 0.28)
  wrong_pi_sum <- c(0.5, 0.7)
  wrong_pi_value <- c(-0.5)

  #  Wrong alpha
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = wrong_alpha_bernoulli,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = wrong_alpha_negative,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = wrong_alpha_negative,
    distribution = "poisson"
  ))
  # No error
  expect_no_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_no_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = alpha,
    distribution = "poisson"
  ))
  expect_no_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi_sum_float, rho = pi_sum_float, alpha = alpha_3,
    distribution = "bernoulli"
  ))

  #  Wrong pi
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = wrong_pi_sum, rho = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = wrong_pi_value, rho = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  #  Wrong rho
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = wrong_pi_sum, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = wrong_pi_value, alpha = alpha,
    distribution = "bernoulli"
  ))

  # No error
  expect_no_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_no_error(colSBM:::generate_bipartite_network(
    nr = nr, nc = nc,
    pi = pi, rho = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
})

test_that("Testing various wrong arguments for generating unipartite networks", {
  n <- 10L

  alpha <- matrix(c(1L, 1L, 1L, 1L), nrow = 2)
  alpha_3 <- matrix(c(
    0.9, 0.4, 0.1,
    0.4, 0.1, 0.05,
    0.1, 0.05, 0.05
  ), byrow = TRUE, ncol = 3)
  wrong_alpha_bernoulli <- matrix(2L)
  wrong_alpha_negative <- matrix(-5L)

  pi <- c(0.5, 0.5)
  pi_sum_float <- c(1 / 3 + 0.28, 1 / 3, 1 / 3 - 0.28)
  wrong_pi_sum <- c(0.5, 0.7)
  wrong_pi_value <- c(-0.5)

  #  Wrong alpha
  expect_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = wrong_alpha_bernoulli,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = wrong_alpha_negative,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = wrong_alpha_negative,
    distribution = "poisson"
  ))

  # No error
  expect_no_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_no_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = alpha,
    distribution = "poisson"
  ))
  expect_no_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi_sum_float, alpha = alpha_3,
    distribution = "bernoulli"
  ))

  #  Wrong pi
  expect_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = wrong_pi_sum, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = wrong_pi_value, alpha = alpha,
    distribution = "bernoulli"
  ))

  # No error
  expect_no_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
  expect_no_error(colSBM:::generate_unipartite_network(
    n = n,
    pi = pi, alpha = alpha,
    distribution = "bernoulli"
  ))
})

test_that("Testing that all arguments work for generating unipartite networks", {
  isNested <- function(l) {
    stopifnot(is.list(l))
    for (i in l) {
      if (is.list(i)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  n <- 10L
  alpha <- matrix(c(1L, 1L, 1L, 1L), nrow = 2)
  pi <- c(0.5, 0.5)

  # Wrong M error
  expect_error(colSBM::generate_unipartite_collection(
    n = c(50, 50, 40),
    pi = pi, alpha = alpha,
    distribution = "bernoulli",
    M = 2
  ))


  expect_false(
    isNested(
      colSBM::generate_unipartite_collection(
        n = n,
        pi = pi,
        alpha = alpha,
        M = 3L
      )
    )
  )

  expect_true(
    isNested(
      colSBM::generate_unipartite_collection(
        n = n,
        pi = pi,
        alpha = alpha,
        M = 3L,
        return_memberships = TRUE
      )
    )
  )

  expect_no_error(colSBM::generate_unipartite_collection(
    n = c(50, 50, 40),
    pi = pi, alpha = alpha,
    distribution = "bernoulli",
    M = 3L
  ))
})

test_that("Base case spectral biclustering", {
  X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
  expect_equal(
    colSBM:::spectral_biclustering(X = X, K = c(1, 1)),
    list(
      row_clustering = c(1, 1),
      col_clustering = c(1, 1)
    )
  )
})

test_that("Spectral clustering too many clusters works", {
  X <- colSBM:::generate_unipartite_network(n = 10, pi = 1, alpha = 0.9)
  expect_message(
    colSBM:::spectral_clustering(X = X, K = 12),
  )
})

test_that("Base case spectral biclustering", {
  X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
  expect_equal(
    colSBM:::spectral_biclustering(X = X, K = c(1, 1)),
    list(
      row_clustering = c(1, 1),
      col_clustering = c(1, 1)
    )
  )
})

test_that("Base case hierarchical biclustering", {
  X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
  expect_equal(
    colSBM:::bipartite_hierarchic_clustering(X = X, K = c(1, 1)),
    list(
      row_clustering = c(1, 1),
      col_clustering = c(1, 1)
    )
  )
})

test_that("Hierarchical biclustering works", {
  X <- colSBM:::generate_bipartite_network(
    nr = 10, 10, pi = c(0.4, 0.6), rho = 1,
    alpha = matrix(c(0.9, 0.1), byrow = TRUE, nrow = 2)
  )
  expect_no_error(
    colSBM:::bipartite_hierarchic_clustering(X = X, K = c(2, 1)),
  )
})

test_that("split_clust fails with a table", {
  X <- as.table(matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1, 2, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 1, 0, 1, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0
  ), nrow = 10, ncol = 10, byrow = TRUE))
  Z <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  expect_no_error(colSBM:::split_clust(X = X, Q = 2, Z = Z, is_bipartite = TRUE))
})

test_that("Helper functions work", {
  expect_identical(colSBM:::logistic(Inf), 1)
  expect_identical(colSBM:::logistic(-Inf), 0)

  expect_identical(colSBM:::logit(0), -Inf)
  expect_identical(colSBM:::logit(1), Inf)
})

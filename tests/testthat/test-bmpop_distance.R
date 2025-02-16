pi1 <- c(0.3, 0.6, 0.1)
pi2 <- c(0.2, 0.3, 0.5)

rho1 <- c(0.1, 0.8, 0.1)
rho2 <- c(1L / 3L, 1L / 3L, 1L / 3L)

alpha1 <- matrix(
  c(
    0.9, 0.4, 0.05,
    0.45, 0.25, 0.05,
    0.1, 0.05, 0.05
  ),
  nrow = 3L
)

alpha2 <- matrix(
  c(
    0.9, 0.05, 0.05,
    0.05, 0.7, 0.05,
    0.05, 0.05, 0.65
  ),
  nrow = 3L
)
# colSBM
test_that("SBM : Distance for same parameters is zero", {
  expect_identical(
    dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
      directed = TRUE
    ), 0
  )
  expect_identical(
    dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L1",
      directed = TRUE
    ), 0
  )
  expect_identical(
    dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L2",
      directed = TRUE
    ), 0
  )
  expect_identical(
    dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L1",
      directed = TRUE
    ), 0
  )
})

test_that("SBM : If directed is missing can compute it", {
  # Undirected
  expect_no_error(
    dist_bmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    dist_bmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  # Mixed
  # Directed
  expect_no_error(
    dist_bmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
})

test_that("SBM : If pi is missing can compute anyway only with alpha", {
  # Undirected
  expect_no_error(
    dist_bmpop_max(
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    dist_bmpop_max(
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  # Mixed
  # Directed
  expect_no_error(
    dist_bmpop_max(
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
})

# colBiSBM
test_that("BiSBM : Distance for same parameters is zero", {
  expect_identical(
    dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    ), 0
  )
  expect_identical(
    dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L1",
    ), 0
  )
  expect_identical(
    dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L2",
    ), 0
  )
  expect_identical(
    dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L1",
    ), 0
  )
})

test_that("BiSBM : If pi or rho is missing can compute anyway only with alpha", {
  # Undirected
  expect_no_error(
    dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  # Mixed
  # Directed
  expect_no_error(
    dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
})

# Graphon distance
## New parameters to account for LBM structure
pi1 <- c(0.3, 0.6, 0.1)
pi2 <- c(0.2, 0.3, 0.5)

# Unipartite compatible
alpha1 <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.7, 0.05,
  0.05, 0.05, 0.65
), nrow = 3L)
alpha2 <- matrix(c(
  0.9, 0.4, 0.05,
  0.45, 0.25, 0.05,
  0.1, 0.05, 0.05
), nrow = 3L)

# Bipartite only
rho1 <- c(0.1, 0.3, 0.2, 0.4)
rho2 <- c(0.2, 0.3, 0.5)

alpha_bip1 <- matrix(c(
  0.9, 0.05, 0.05, 0.3,
  0.05, 0.7, 0.05, 0.05,
  0.05, 0.05, 0.65, 0.05
), nrow = 3L)
alpha_bip2 <- matrix(c(
  0.9, 0.4, 0.05,
  0.45, 0.25, 0.05,
  0.1, 0.05, 0.05
), nrow = 3L)

# Tests for graphon_distance_bipartite
test_that("graphon_distance_bipartite computes correct distance", {
  result <- graphon_distance_bipartite(list(pi1, pi2), list(rho1, rho2), list(alpha_bip1, alpha_bip2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for dist_graphon_bipartite_all_permutations
test_that("dist_graphon_bipartite_all_permutations computes correct distance", {
  result <- dist_graphon_bipartite_all_permutations(list(pi1, pi2), list(rho1, rho2), list(alpha_bip1, alpha_bip2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for dist_graphon_bipartite_marginals
test_that("dist_graphon_bipartite_marginals computes correct distance", {
  result <- dist_graphon_bipartite_marginals(list(pi1, pi2), list(rho1, rho2), list(alpha_bip1, alpha_bip2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for dist_graphon_bipartite_symmetrization
test_that("dist_graphon_bipartite_symmetrization computes correct distance", {
  result <- dist_graphon_bipartite_symmetrization(list(pi1, pi2), list(rho1, rho2), list(alpha_bip1, alpha_bip2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for graphon_distance_unipartite
test_that("graphon_distance_unipartite computes correct distance", {
  result <- graphon_distance_unipartite(list(pi1, pi2), list(alpha1, alpha2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for dist_graphon_unipartite_all_permutations
test_that("dist_graphon_unipartite_all_permutations computes correct distance", {
  result <- dist_graphon_unipartite_all_permutations(list(pi1, pi2), list(alpha1, alpha2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

# Tests for dist_graphon_unipartite_marginals
test_that("dist_graphon_unipartite_marginals computes correct distance", {
  result <- dist_graphon_unipartite_marginals(list(pi1, pi2), list(alpha1, alpha2))
  expect_type(result, "double")
  expect_gte(result, 0)
})

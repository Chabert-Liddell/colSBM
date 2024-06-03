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
    colSBM::dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
      directed = TRUE
    ), 0
  )
  expect_identical(
    colSBM::dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L1",
      directed = TRUE
    ), 0
  )
  expect_identical(
    colSBM::dist_bmpop_max(
      pi = list(pi1, pi1),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L2",
      directed = TRUE
    ), 0
  )
  expect_identical(
    colSBM::dist_bmpop_max(
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
    colSBM::dist_bmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    colSBM::dist_bmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  #  Mixed
  # Directed
  expect_no_error(
    colSBM::dist_bmpop_max(
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
    colSBM::dist_bmpop_max(
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    colSBM::dist_bmpop_max(
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  #  Mixed
  # Directed
  expect_no_error(
    colSBM::dist_bmpop_max(
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
})

#  colBiSBM
test_that("BiSBM : Distance for same parameters is zero", {
  expect_identical(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    ), 0
  )
  expect_identical(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L1",
    ), 0
  )
  expect_identical(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi1),
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "mean",
      norm = "L2",
    ), 0
  )
  expect_identical(
    colSBM::dist_bisbmpop_max(
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
    colSBM::dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha1),
      weight = "max",
      norm = "L2",
    )
  )
  # Directed
  expect_no_error(
    colSBM::dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha2, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  #  Mixed
  # Directed
  expect_no_error(
    colSBM::dist_bisbmpop_max(
      rho = list(rho1, rho2),
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
  expect_no_error(
    colSBM::dist_bisbmpop_max(
      pi = list(pi1, pi2),
      alpha = list(alpha1, alpha2),
      weight = "max",
      norm = "L2",
    )
  )
})

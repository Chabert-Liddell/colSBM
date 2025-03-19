colsbm <- colSBM::bmpop$new(netlist = foodwebs[1:2], directed = TRUE, distribution = "poisson")
test_that("lfactorial", {
  expect_equal(colsbm$logfactA[1], 0)
})

test_that("split_clust", {
  expect_equal(
    sort(
      unique(
        colSBM:::split_clust(
          as.matrix(colsbm$A[[1]]), rep(1, nrow(colsbm$A[[1]])), 1
        )[[1]]
      )
    ),
    c(1, 2)
  )
})

test_that("Spectral_init for missing models works", {
  cloned_colsbm <- colsbm$clone()
  cloned_colsbm$optimize()
  cloned_colsbm$model_list[[1]][[2]] <- NULL
  expect_no_error(cloned_colsbm$forward_pass(Q_min = 1, Q_max = 4))
})

test_that("clusterize_unipartite_networks_iid", {
  set.seed(1234)
  Net <- lapply(
    list(.7, .7, .2, .2),
    function(p) {
      A <- matrix(0, 15, 15)
      A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
      A <- A + t(A)
    }
  )
  cl <- clusterize_unipartite_networks(Net,
    colsbm_model = "iid",
    directed = FALSE,
    distribution = "bernoulli",
    nb_run = 1,
    global_opts = list(nb_cores = 1),
    fit_opts = list(
      Q_max = 2,
      nb_init = 1,
      nb_models = 1,
      depth = 1
    )
  )
  expect_identical(
    {
      cl$partition[[1]]$M
    },
    4L
  )
})

test_that("estimate_colsbm_poisson_delta with verbosity", {
  expect_equal(
    {
      set.seed(1234)
      Net <- lapply(
        list(.7, .7, .2, .2),
        function(p) {
          A <- matrix(0, 15, 15)
          A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
          A <- A + t(A)
        }
      )
      fit <- estimate_colSBM(Net,
        colsbm_model = "delta",
        directed = FALSE,
        distribution = "poisson",
        nb_run = 1,
        global_opts = list(
          nb_cores = 1L,
          verbosity = 4L
        ),
        fit_opts = list(
          Q_max = 2,
          nb_init = 1,
          nb_models = 1,
          depth = 1
        )
      )
      fit$best_fit$Q
    },
    1
  )
  expect_no_error(fit$show())
})

test_that("estimate_colsbm_poisson_delta NULL fit_sbm", {
  expect_equal(
    {
      set.seed(1234)
      Net <- lapply(
        list(.7, .7, .2, .2),
        function(p) {
          A <- matrix(0, 15, 15)
          A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
          A <- A + t(A)
        }
      )
      fit <- estimate_colSBM(Net,
        colsbm_model = "delta",
        directed = FALSE,
        distribution = "poisson",
        nb_run = 1,
        fit_sbm = NULL,
        global_opts = list(nb_cores = 1),
        fit_opts = list(
          Q_max = 2,
          nb_init = 1,
          nb_models = 1,
          depth = 1
        )
      )
      fit$best_fit$Q
    },
    1
  )
})

test_that("estimate_colsbm_poisson_delta with NA", {
  expect_equal(
    {
      set.seed(1234)
      Net <- lapply(
        list(.7, .7, .2, .2),
        function(p) {
          A <- matrix(0, 15, 15)
          A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
          A <- A + t(A)
          A[sample(x = nrow(A), size = nrow(A) / 5), sample(x = ncol(A), size = ncol(A) / 5L)] <- NA
          A
        }
      )
      fit <- estimate_colSBM(Net,
        colsbm_model = "delta",
        directed = FALSE,
        distribution = "poisson",
        nb_run = 1,
        global_opts = list(nb_cores = 1),
        fit_opts = list(
          Q_max = 2,
          nb_init = 1,
          nb_models = 1,
          depth = 1
        )
      )
      fit$best_fit$Q
    },
    1
  )
})

test_that("estimate_colsbm_poisson with pi model", {
  expect_equal(
    {
      set.seed(1234)
      Net <- lapply(
        list(.7, .7, .7, .7),
        function(p) {
          A <- matrix(0, 15, 15)
          A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
          A <- A + t(A)
          A
        }
      )
      fit <- estimate_colSBM(Net,
        colsbm_model = "pi",
        directed = FALSE,
        distribution = "poisson",
        nb_run = 1,
        global_opts = list(nb_cores = 1),
        fit_opts = list(
          Q_max = 2,
          nb_init = 1,
          nb_models = 1,
          depth = 1
        )
      )
      fit$best_fit$Q
    },
    1
  )
})

test_that("estimate_colsbm_poisson with deltapi model", {
  expect_equal(
    {
      set.seed(1234)
      Net <- lapply(
        list(.7, .7, .7, .7),
        function(p) {
          A <- matrix(0L, 15L, 15L)
          A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1L
          A <- A + t(A)
          A
        }
      )
      fit <- estimate_colSBM(Net,
        colsbm_model = "deltapi",
        directed = FALSE,
        distribution = "poisson",
        nb_run = 1L,
        global_opts = list(nb_cores = 1L),
        fit_opts = list(
          Q_max = 2L,
          nb_init = 1L,
          nb_models = 1L,
          depth = 1L
        )
      )
      fit$best_fit$Q
    },
    1
  )
})

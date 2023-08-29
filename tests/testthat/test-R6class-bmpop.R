colsbm <- colSBM::bmpop$new(netlist = foodwebs[1:2], directed = TRUE, distribution = "poisson")
test_that("lfactorial", {
  expect_equal(colsbm$logfactA[1], 0)
})
test_that("split_clust", {
  expect_equal(
    sort(unique(split_clust(as.matrix(colsbm$A[[1]]), rep(1, nrow(colsbm$A[[1]])), 1)[[1]])),
    c(1, 2)
  )
})

test_that("clusterize_networks_iid", {
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
      cl <- clusterize_networks(Net,
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
      c(diff(cl[[2]]$net_id), diff(cl[[3]]$net_id))
    },
    c(1, 1)
  )
})
test_that("estimate_colsbm_poisson_delta", {
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

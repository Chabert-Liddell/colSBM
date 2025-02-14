# colSBM
test_that("colsbm-model not correct raises an error", {
  set.seed(1234L)
  expect_error(estimate_colSBM(list(matrix(runif(100, 10, 10))),
    netmodel = "pid"
  ))
})

test_that("Estimate colSBM runs without problems", {
  set.seed(1234L)
  expect_no_error(estimate_colSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model =
      "iid",
    global_opts = list(
      nb_cores = 2L,
      backend = "no_mc",
      verbosity = 0L
    )
  ))
})

test_that("Wrong model raises an error", {
  set.seed(1234L)
  expect_error(estimate_colSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model =
      "wrongmodel",
    global_opts = list(
      nb_cores = 2L,
      backend = "no_mc",
      verbosity = 0L
    )
  ))
})

test_that("Estimate network benefitting from joint modelisation with null nb_cores", {
  set.seed(1234L)
  iid_collection <- generate_unipartite_collection(
    n = 50,
    pi = c(0.4, 0.6),
    alpha = matrix(
      c(
        0.9, 0.05,
        0.05, 0.4
      ),
      nrow = 2L
    ),
    M = 3
  )
  fit <- estimate_colSBM(
    netlist = iid_collection,
    colsbm_model = "iid",
    nb_run = 1L,
    global_opts = list(verbosity = 4L, nb_cores = NULL, Q_max = 5L, Q_min = 1L)
  )
  expect_identical(fit$best_fit$Q, 2L)
})

# colBiSBM
test_that("colbisbm-model not correct raises an error", {
  set.seed(1234L)
  expect_error(estimate_colBiSBM(
    netlist = list(matrix(runif(100, 10, 10))),
    colsbm_model = "pirhod"
  ))
})

test_that("Estimate colBiSBM does not accept incorrect nb_cores", {
  set.seed(1234L)
  expect_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 0L,
      verbosity = 0L
    )
  ))
})

test_that("Estimate colBiSBM does not authorize non binary data for bernoulli distribution", {
  set.seed(1234L)
  expect_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(rpois(100, lambda = 5)), nrow = 10),
      matrix(as.numeric(rpois(100, lambda = 5)), nrow = 10)
    ),
    distribution = "bernoulli",
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 1L,
      verbosity = 0L
    )
  ))
})

test_that("Estimate colBiSBM does not authorize non integer data for poisson distribution", {
  set.seed(1234L)
  expect_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100)), nrow = 10),
      matrix(as.numeric(runif(100)), nrow = 10)
    ),
    distribution = "poisson",
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 1L,
      verbosity = 0L
    )
  ))
})

test_that("Estimate colBiSBM does not authorize negative integers for poisson distribution", {
  set.seed(1234L)
  expect_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(-rpois(100, 5)), nrow = 10),
      matrix(as.numeric(-rpois(100, 5)), nrow = 10)
    ),
    distribution = "poisson",
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 1L,
      verbosity = 0L
    )
  ))
})

test_that("Estimate colBiSBM runs without problems with verbosity", {
  set.seed(1234L)
  expect_no_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 2L,
      backend = "no_mc",
      verbosity = 4L,
      Q1_max = NULL,
      Q2_max = NULL
    )
  ))
})

test_that("Estimate colBiSBM runs with pi, rho and pirho. With and without NA", {
  set.seed(1234L)
  pirho_collection <- generate_bipartite_collection(
    nr = 50,
    nc = 50,
    pi = c(0.1, 0.3, 0.6),
    rho = c(0.8, 0.2),
    alpha = matrix(
      c(
        0.2, 0.9,
        0.05, 0.5,
        0.05, 0.05
      ),
      nrow = 3L
    ),
    M = 3L,
    model = "pirho"
  )
  pirho_collection_na <- pirho_collection
  NA_index <- sample.int(length(pirho_collection_na[[1]]), size = 20)
  pirho_collection_na[[1]][NA_index] <- NA
  expect_no_error(
    fit_iid <- estimate_colBiSBM(
      netlist = pirho_collection_na,
      colsbm_model = "iid",
      nb_run = 2L,
      global_opts = list(
        nb_cores = 2L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_iid <- estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "iid",
      nb_run = 2L,
      global_opts = list(
        nb_cores = 2L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_pi <- estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "pi",
      global_opts = list(
        nb_cores = 2L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_rho <- estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "rho",
      global_opts = list(
        nb_cores = 2L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_pirho <- estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "pirho",
      global_opts = list(
        nb_cores = 2L,
        backend = "no_mc",
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )

  expect_warning(adjust_colBiSBM(
    fitted_bisbmpop = fit_iid,
    Q = c(1L, 1L)
  ))
})


test_that("Check that LBM and colBiSBM with one network return the same", {
  set.seed(1234)
  testNet <- generate_bipartite_network(
    nr = 60L, nc = 60L,
    pi = c(0.3, 0.7),
    rho = c(0.45, 0.55),
    alpha = matrix(c(
      0.8, 0.4,
      0.2, 0.05
    ), nrow = 2L)
  )
  init_lbm <- sbm::estimateBipartiteSBM(
    netMat = testNet,
    model = "bernoulli",
    estimOptions = list(
      verbosity = 0L,
      plot = 0L,
      nbcores = 1L
    )
  )
  lbm <- fitBipartiteSBMPop[["new"]](
    A = list(testNet),
    Q = init_lbm[["nbBlocks"]],
    distribution = "bernoulli",
    free_mixture_row = FALSE,
    free_mixture_col = FALSE,
    init_method = "given",
    Z = list(init_lbm[["memberships"]]),
    fit_opts = list(verbosity = 0L)
  )
  lbm[["optimize"]]()

  colBiSBM <- estimate_colBiSBM(
    netlist = list(testNet), colsbm_model = "iid",
    nb_run = 1L,
    global_opts = list(verbosity = 0L, backend = "no_mc")
  )

  collbm <- colBiSBM[["best_fit"]]

  expect_equal(unname(lbm[["ICL"]]), collbm[["ICL"]], tolerance = 1e-3)
  expect_equal(unname(lbm[["entropy"]]), collbm[["entropy"]], tolerance = 1e-3)
  expect_identical(unname(lbm[["penalty"]]), collbm[["penalty"]])
  expect_equal(unname(lbm[["BICL"]]), collbm[["BICL"]], tolerance = 1e-3)
})

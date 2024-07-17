# colSBM
test_that("colsbm-model not correct raises an error", {
  set.seed(1234L)
  expect_error(colSBM::estimate_colSBM(list(matrix(runif(100, 10, 10))),
    netmodel = "pid"
  ))
})

test_that("Estimate colSBM runs without problems", {
  set.seed(1234L)
  expect_no_error(colSBM::estimate_colSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model =
      "iid",
    global_opts = list(
      nb_cores = 2L,
      verbosity = 0L
    )
  ))
})

test_that("Wrong model raises an error", {
  set.seed(1234L)
  expect_error(colSBM::estimate_colSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model =
      "wrongmodel",
    global_opts = list(
      nb_cores = 2L,
      verbosity = 0L
    )
  ))
})

test_that("Estimate network benefitting from joint modelisation with null nb_cores", {
  set.seed(1234L)
  iid_collection <- colSBM::generate_unipartite_collection(
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
  fit <- colSBM::estimate_colSBM(
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
  expect_error(colSBM::estimate_colBiSBM(
    netlist = list(matrix(runif(100, 10, 10))),
    colsbm_model = "pirhod"
  ))
})

test_that("Estimate colBiSBM does not accept incorrect nb_cores", {
  set.seed(1234L)
  expect_error(colSBM::estimate_colBiSBM(
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
  expect_error(colSBM::estimate_colBiSBM(
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

test_that("Estimate colBiSBM runs without problems with verbosity", {
  set.seed(1234L)
  expect_no_error(colSBM::estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model = "iid",
    global_opts = list(
      nb_cores = 2L,
      verbosity = 4L,
      Q1_max = NULL,
      Q2_max = NULL
    )
  ))
})

test_that("Estimate colBiSBM runs with pi, rho and pirho. With and without NA", {
  set.seed(1234L)
  pirho_collection <- colSBM::generate_bipartite_collection(
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
    fit_iid <- colSBM::estimate_colBiSBM(
      netlist = pirho_collection_na,
      colsbm_model = "iid",
      nb_run = 10L,
      global_opts = list(
        nb_cores = 2L,
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_iid <- colSBM::estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "iid",
      nb_run = 10L,
      global_opts = list(
        nb_cores = 2L,
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_pi <- colSBM::estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "pi",
      global_opts = list(
        nb_cores = 2L,
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_rho <- colSBM::estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "rho",
      global_opts = list(
        nb_cores = 2L,
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )
  expect_no_error(
    fit_pirho <- colSBM::estimate_colBiSBM(
      netlist = pirho_collection,
      colsbm_model = "pirho",
      global_opts = list(
        nb_cores = 2L,
        verbosity = 0L,
        Q1_max = 10L,
        Q2_max = 10L
      )
    )
  )

  expect_warning(colSBM::adjust_colBiSBM(
    fitted_bisbmpop = fit_iid,
    Q = c(1L, 1L)
  ))
})

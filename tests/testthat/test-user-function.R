test_that("colsbm-model not correct raises an error", {
  expect_error(estimate_colSBM(list(matrix(runif(100, 10, 10))),
    netmodel = "pid"
  ))
})

test_that("colbisbm-model not correct raises an error", {
  expect_error(estimate_colBiSBM(
    netlist = list(matrix(runif(100, 10, 10))),
    colsbm_model = "pirhod"
  ))
})

test_that("Estimate colBiSBM does not accept incorrect nb_cores", {
  expect_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model = "iid",
    global_opts = list(nb_cores = 0L)
  ))
})

test_that("Estimate colBiSBM runs without problems", {
  expect_no_error(estimate_colBiSBM(
    netlist = list(
      matrix(as.numeric(runif(100) > 0.5), nrow = 10),
      matrix(as.numeric(runif(100) > 0.5), nrow = 10)
    ),
    colsbm_model = "iid",
    global_opts = list(nb_cores = 2L)
  ))
})
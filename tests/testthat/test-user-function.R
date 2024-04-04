test_that("colsbm-model not correct raises an error", {
  expect_error(estimate_colSBM(list(matrix(runif(100, 10, 10))),
    netmodel = "pid"
  ))
})

test_that("colbisbm-model not correct raises an error", {
  expect_error(estimate_colBiSBM(netlist = list(matrix(runif(100, 10, 10))),
    colsbm_model = "pirhod"
  ))
})
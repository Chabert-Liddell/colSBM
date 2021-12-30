test_that("colsbm-model", {
  expect_error(estimate_colSBM(list(matrix(runif(100, 10, 10))),
                               netmodel  = "pid"))
})

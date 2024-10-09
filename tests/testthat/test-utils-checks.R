# check_is_integer_over_thresh
test_that("check_is_integer_over_thresh() detects incorrect type input", {
  expect_error(check_is_integer_over_thresh("a", 1L), "must be an integer")
  expect_error(check_is_integer_over_thresh(1.5, 1L), "must be an integer")
  expect_error(check_is_integer_over_thresh(matrix(1, 2), 1L), "must be an integer")
})

test_that("check_is_integer_over_thresh() detects input below threshold", {
  expect_error(check_is_integer_over_thresh(0L, 1L), "must be at least 1")
  expect_error(check_is_integer_over_thresh(0L, 10L), "must be at least 10")
})

test_that("check_is_integer_over_thresh() works correctly", {
  expect_identical(check_is_integer_over_thresh(1L, 1L), 1L)
  expect_identical(check_is_integer_over_thresh(10L, 1L), 10L)
})

# check_networks_list
test_that("check_networks_list() detects empty input", {
  expect_error(check_networks_list(arg = NULL), "must be supplied")
})
test_that("check_networks_list() detects incorrect type input", {
  expect_error(check_networks_list("a"), "must be a list of matrices")
  expect_error(check_networks_list(list()), "must be a list of matrices")
  expect_error(check_networks_list(list(matrix(1, 2), 1)), "must be a list of matrices")
})

test_that("check_networks_list() detects input below threshold", {
  expect_error(check_networks_list(list(matrix(1, 2)), 2L), "must be of length at least 2")
  expect_error(check_networks_list(list(matrix(1, 2), matrix(1, 2)), 3L), "must be of length at least 3")
})

test_that("check_networks_list() works correctly", {
  expect_no_error(check_networks_list(list(matrix(1, 2), matrix(1, 2)), 2L))
  expect_no_error(check_networks_list(list(matrix(1, 2), matrix(1, 2), matrix(1, 2)), 2L))
})

# check_dissimilarity_matrix
test_that("check_dissimilarity_matrix() detects empty input", {
  expect_error(check_dissimilarity_matrix(arg = NULL), "must be supplied")
})

test_that("check_dissimilarity_matrix() detects incorrect type input", {
  expect_error(check_dissimilarity_matrix("a"), "must be a dissimilarity matrix")
  expect_error(check_dissimilarity_matrix(matrix(c("a", "b"), 2)), "must be a dissimilarity matrix")
  expect_error(check_dissimilarity_matrix(matrix(c(TRUE, FALSE))), "must be a dissimilarity matrix")
})

test_that("check_dissimilarity_matrix() works correctly", {
  expect_no_error(check_dissimilarity_matrix(matrix(1, 2)))
  expect_no_error(check_dissimilarity_matrix(matrix(c(1, 2, 3, 4), 2)))
  expect_no_error(check_dissimilarity_matrix(matrix(c(1, 2, 3, 4), 2, dimnames = list(letters[1:2], letters[1:2]))))
})

# check_bipartite_colsbm_models
test_that("check_bipartite_colsbm_models() detects empty input", {
  expect_error(check_bipartite_colsbm_models(arg = NULL), "must be supplied")
})

test_that("check_bipartite_colsbm_models() detects incorrect choice", {
  expect_error(check_bipartite_colsbm_models("a"), "must be one of")
  expect_error(check_bipartite_colsbm_models(1), "must be a character vector")
  expect_error(check_bipartite_colsbm_models(c("a", "b")), "must be one of")
})

test_that("check_bipartite_colsbm_models() works correctly", {
  expect_no_error(check_bipartite_colsbm_models("iid"))
  expect_no_error(check_bipartite_colsbm_models("pi"))
  expect_no_error(check_bipartite_colsbm_models("rho"))
  expect_no_error(check_bipartite_colsbm_models("pirho"))
})

# check_unipartite_colsbm_models
test_that("check_unipartite_colsbm_models() detects empty input", {
  expect_error(check_unipartite_colsbm_models(arg = NULL), "must be supplied")
})

test_that("check_unipartite_colsbm_models() detects incorrect choice", {
  expect_error(check_unipartite_colsbm_models("a"), "must be one of")
  expect_error(check_unipartite_colsbm_models(1), "must be a character vector")
  expect_error(check_unipartite_colsbm_models(c("a", "b")), "must be one of")
})

test_that("check_unipartite_colsbm_models() works correctly", {
  expect_no_error(check_unipartite_colsbm_models("iid"))
  expect_no_error(check_unipartite_colsbm_models("pi"))
  expect_no_error(check_unipartite_colsbm_models("delta"))
  expect_no_error(check_unipartite_colsbm_models("deltapi"))
})

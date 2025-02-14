common_fit_opts <- list(max_vem_steps = 3000L)

test_that("clusterize_bipartite_networks works with valid inputs", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    global_opts = list(nb_cores = 1),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result, inherits, "fitBipartiteSBMPop")))
})

test_that("clusterize_bipartite_networks handles invalid colsbm_model", {
  netlist <- list(matrix(0, 10, 10), matrix(0, 10, 10))

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "invalid_model",
      global_opts = list(nb_cores = 1),
      fit_opts = common_fit_opts
    ),
    "`colsbm_model` must be one of"
  )
})

test_that("clusterize_bipartite_networks handles empty netlist", {
  netlist <- list()

  expect_error(
    clusterize_bipartite_networks(
      netlist = netlist,
      colsbm_model = "iid",
      global_opts = list(nb_cores = 1),
      fit_opts = common_fit_opts
    ),
    "`netlist` must be a list of matrices."
  )
})

test_that("clusterize_bipartite_networks works with different distributions", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    distribution = "poisson",
    global_opts = list(nb_cores = 1),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result, inherits, "fitBipartiteSBMPop")))
})

test_that("clusterize_bipartite_networks works with full_inference = TRUE", {
  set.seed(0)
  alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
  alpha2 <- matrix(c(0.8, 0.5, 0.5, 0.2), byrow = TRUE, nrow = 2)
  first_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha1, M = 2
  )
  second_collection <- generate_bipartite_collection(
    nr = 50, nc = 25,
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    alpha = alpha2, M = 2
  )
  netlist <- append(first_collection, second_collection)

  result <- clusterize_bipartite_networks(
    netlist = netlist,
    colsbm_model = "iid",
    full_inference = TRUE,
    global_opts = list(nb_cores = 1),
    fit_opts = common_fit_opts
  )

  expect_type(result, "list")
  expect_true(all(sapply(result, inherits, "fitBipartiteSBMPop")))
})

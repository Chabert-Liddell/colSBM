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

# check_colsbm_emission_distribution
test_that("check_colsbm_emission_distribution() detects empty input", {
  expect_error(check_colsbm_emission_distribution(arg = NULL), "must be supplied")
})

test_that("check_colsbm_emission_distribution() detects incorrect choice", {
  expect_error(check_colsbm_emission_distribution("a"), "must be one of")
  expect_error(check_colsbm_emission_distribution(1), "must be a character vector")
  expect_error(check_colsbm_emission_distribution(c("a", "b")), "must be one of")
})

test_that("check_colsbm_emission_distribution() works correctly", {
  expect_no_error(check_colsbm_emission_distribution("bernoulli"))
  expect_no_error(check_colsbm_emission_distribution("poisson"))
})

# check_networks_list_match_emission_distribution
test_that("check_networks_list_match_emission_distribution() detects empty input", {
  expect_error(check_networks_list_match_emission_distribution(arg = NULL, distrib_arg = NULL), "must be supplied")
})

test_that("check_networks_list_match_emission_distribution() detects incorrect choice", {
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1))), emission_distribution = "a"), "must be one of")
})

#  Bernoulli
## Non matching networks matrices
test_that("check_networks_list_match_emission_distribution() detects non matching networks matrices for bernoulli", {
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0.5, 1))), emission_distribution = "bernoulli"), "Non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1.5))), emission_distribution = "bernoulli"), "Non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1, "0"))), emission_distribution = "bernoulli"), "Non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1, 2))), emission_distribution = "bernoulli"), "either 0 or 1")
})

## Matching networks matrices
test_that("check_networks_list_match_emission_distribution() detects matching networks matrices for bernoulli", {
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1))), emission_distribution = "bernoulli"))
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1), c(1, 0))), emission_distribution = "bernoulli"))
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1), c(1, 0)), matrix(c(0, 1), c(1, 0))), emission_distribution = "bernoulli"))
})


# Poisson
## Non matching networks matrices
test_that("check_networks_list_match_emission_distribution() detects non matching networks matrices for poisson", {
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0.5, 1))), emission_distribution = "poisson"), "non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1.5))), emission_distribution = "poisson"), "non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1, "0"))), emission_distribution = "poisson"), "non integer entries")
  expect_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, -1, 2))), emission_distribution = "poisson"), "with negative entries")
})

## Matching networks matrices
test_that("check_networks_list_match_emission_distribution() detects matching networks matrices for poisson", {
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 5))), emission_distribution = "poisson"))
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 10), c(2, 3))), emission_distribution = "poisson"))
  expect_no_error(check_networks_list_match_emission_distribution(networks_list = list(matrix(c(0, 1), c(1, 5)), matrix(c(0, 1), c(1, 0))), emission_distribution = "poisson"))
})

# check_net_id
test_that("check_net_id() detects incorrect type input", {
  expect_error(check_net_id(1, list(matrix(1, 2))), "must be a character vector")
  expect_error(check_net_id(TRUE, list(matrix(1, 2))), "must be a character vector")
})

test_that("check_net_id() detects length mismatch", {
  expect_error(check_net_id(c("id1"), list(matrix(1, 2), matrix(1, 2))), "must have the same length as")
  expect_error(check_net_id(c("id1", "id2", "id3"), list(matrix(1, 2), matrix(1, 2))), "must have the same length as")
})

test_that("check_net_id() works correctly", {
  expect_no_error(check_net_id(c("id1", "id2"), list(matrix(1, 2), matrix(1, 2))))
  expect_no_error(check_net_id(c("id1", "id2", "id3"), list(matrix(1, 2), matrix(1, 2), matrix(1, 2))))
})

# Initializing default global options
default_go <- default_global_opts_bipartite(netlist = list(matrix(1)))

# check_global_opts
test_that("check_global_opts() detects incorrect type input", {
  expect_error(check_global_opts("a"), "must be a list")
  expect_error(check_global_opts(1), "must be a list")
  expect_error(check_global_opts(TRUE), "must be a list")
})

go1 <- default_go
go1$nb_cores <- "a"
test_that("check_global_opts() detects incorrect nb_cores", {
  expect_error(check_global_opts(go1), "must be an integer")
})

go2 <- default_go
go2$nb_cores <- 1.5
test_that("check_global_opts() detects incorrect nb_cores", {
  expect_error(check_global_opts(go2), "must be an integer")
})

go3 <- default_go
go3$nb_cores <- 0L
test_that("check_global_opts() detects incorrect nb_cores", {
  expect_error(check_global_opts(go3), "must be at least 1")
})

go4 <- default_go
go4$Q1_max <- "a"
test_that("check_global_opts() detects incorrect Q1_max", {
  expect_error(check_global_opts(go4), "must be an integer")
})

go5 <- default_go
go5$Q1_max <- 1.5
test_that("check_global_opts() detects incorrect Q1_max", {
  expect_error(check_global_opts(go5), "must be an integer")
})

go6 <- default_go
go6$Q1_max <- 0L
test_that("check_global_opts() detects incorrect Q1_max", {
  expect_error(check_global_opts(go6), "must be at least 1")
})

go7 <- default_go
go7$Q2_max <- "a"
test_that("check_global_opts() detects incorrect Q2_max", {
  expect_error(check_global_opts(go7), "must be an integer")
})

go8 <- default_go
go8$Q2_max <- 1.5
test_that("check_global_opts() detects incorrect Q2_max", {
  expect_error(check_global_opts(go8), "must be an integer")
})

go9 <- default_go
go9$Q2_max <- 0L
test_that("check_global_opts() detects incorrect Q2_max", {
  expect_error(check_global_opts(go9), "must be at least 1")
})

go10 <- default_go
go10$nb_init <- "a"
test_that("check_global_opts() detects incorrect nb_init", {
  expect_error(check_global_opts(go10), "must be an integer")
})

go11 <- default_go
go11$nb_init <- 1.5
test_that("check_global_opts() detects incorrect nb_init", {
  expect_error(check_global_opts(go11), "must be an integer")
})

go12 <- default_go
go12$nb_init <- 0L
test_that("check_global_opts() detects incorrect nb_init", {
  expect_error(check_global_opts(go12), "must be at least 1")
})

go13 <- default_go
go13$nb_models <- "a"
test_that("check_global_opts() detects incorrect nb_models", {
  expect_error(check_global_opts(go13), "must be an integer")
})

go14 <- default_go
go14$nb_models <- 1.5
test_that("check_global_opts() detects incorrect nb_models", {
  expect_error(check_global_opts(go14), "must be an integer")
})

go15 <- default_go
go15$nb_models <- 0L
test_that("check_global_opts() detects incorrect nb_models", {
  expect_error(check_global_opts(go15), "must be at least 1")
})

go16 <- default_go
go16$backend <- 1
test_that("check_global_opts() detects incorrect backend", {
  expect_error(check_global_opts(go16), "must be a character vector")
})

go17 <- default_go
go17$depth <- "a"
test_that("check_global_opts() detects incorrect depth", {
  expect_error(check_global_opts(go17), "must be an integer")
})

go18 <- default_go
go18$depth <- 1.5
test_that("check_global_opts() detects incorrect depth", {
  expect_error(check_global_opts(go18), "must be an integer")
})

go19 <- default_go
go19$depth <- 0L
test_that("check_global_opts() detects incorrect depth", {
  expect_error(check_global_opts(go19), "must be at least 1")
})

go20 <- default_go
go20$plot_details <- "a"
test_that("check_global_opts() detects incorrect plot_details", {
  expect_error(check_global_opts(go20), "must be an integer")
})

go21 <- default_go
go21$plot_details <- 1.5
test_that("check_global_opts() detects incorrect plot_details", {
  expect_error(check_global_opts(go21), "must be an integer")
})

go22 <- default_go
go22$plot_details <- -1L
test_that("check_global_opts() detects incorrect plot_details", {
  expect_error(check_global_opts(go22), "must be at least 0")
})

go23 <- default_go
go23$max_pass <- "a"
test_that("check_global_opts() detects incorrect max_pass", {
  expect_error(check_global_opts(go23), "must be an integer")
})

go24 <- default_go
go24$max_pass <- 1.5
test_that("check_global_opts() detects incorrect max_pass", {
  expect_error(check_global_opts(go24), "must be an integer")
})

go25 <- default_go
go25$max_pass <- 0L
test_that("check_global_opts() detects incorrect max_pass", {
  expect_error(check_global_opts(go25), "must be at least 1")
})

go26 <- default_go
go26$verbosity <- "a"
test_that("check_global_opts() detects incorrect verbosity", {
  expect_error(check_global_opts(go26), "must be an integer")
})

go27 <- default_go
go27$verbosity <- 1.5
test_that("check_global_opts() detects incorrect verbosity", {
  expect_error(check_global_opts(go27), "must be an integer")
})

go28 <- default_go
go28$verbosity <- -1L
test_that("check_global_opts() detects incorrect verbosity", {
  expect_error(check_global_opts(go28), "must be at least 0")
})

test_that("check_global_opts() works correctly", {
  expect_no_error(check_global_opts(default_go))
  expect_no_error(check_global_opts(list(nb_cores = 2L, Q1_max = 3L, Q2_max = 4L)))
})

# check_fit_opts
test_that("check_fit_opts() detects incorrect type input", {
  expect_error(check_fit_opts("a"), "must be a list")
  expect_error(check_fit_opts(1), "must be a list")
  expect_error(check_fit_opts(TRUE), "must be a list")
})

fo <- default_fit_opts_bipartite()

fo1 <- fo
fo1$algo_ve <- "invalid"
test_that("check_fit_opts() detects incorrect algo_ve", {
  expect_error(check_fit_opts(fo1), "must be one of")
})

fo2 <- fo
fo2$minibatch <- "not_boolean"
test_that("check_fit_opts() detects incorrect minibatch", {
  expect_error(check_fit_opts(fo2), "must be a boolean")
})

fo3 <- fo
fo3$verbosity <- -1
test_that("check_fit_opts() detects incorrect verbosity", {
  expect_error(check_fit_opts(fo3), "must be at least 0")
})

fo4 <- fo
fo4$tolerance <- "not_double"
test_that("check_fit_opts() detects incorrect tolerance", {
  expect_error(check_fit_opts(fo4), "must be a double")
})

fo5 <- fo
fo5$greedy_exploration_max_steps <- 0
test_that("check_fit_opts() detects incorrect greedy_exploration_max_steps", {
  expect_error(check_fit_opts(fo5), "must be at least 1")
})

fo6 <- fo
fo6$greedy_exploration_max_steps_without_improvement <- 0
test_that("check_fit_opts() detects incorrect greedy_exploration_max_steps_without_improvement", {
  expect_error(check_fit_opts(fo6), "must be at least 1")
})

test_that("check_fit_opts() works correctly", {
  expect_no_error(check_fit_opts(fo))
})

# check_net_id_and_initialize
test_that("check_net_id_and_initialize() initializes net_id correctly", {
  networks_list <- list(matrix(1, 2), matrix(1, 2))
  expect_identical(check_net_id_and_initialize(NULL, networks_list), c(1L, 2L))
})

test_that("check_net_id_and_initialize() detects incorrect type input", {
  networks_list <- list(matrix(1, 2), matrix(1, 2))
  expect_error(check_net_id_and_initialize(TRUE, networks_list), "must be a character or an integer vector")
  expect_error(check_net_id_and_initialize(c(NA, 1), networks_list), "must be a character or an integer vector")
  expect_error(check_net_id_and_initialize(c(NA, "a"), networks_list), "must be a character or an integer vector")
})

test_that("check_net_id_and_initialize() detects length mismatch", {
  networks_list <- list(matrix(1, 2), matrix(1, 2))
  expect_error(check_net_id_and_initialize(c("id1"), networks_list), "must have the same length as")
  expect_error(check_net_id_and_initialize(c("id1", "id2", "id3"), networks_list), "must have the same length as")
})

test_that("check_net_id_and_initialize() works correctly", {
  networks_list <- list(matrix(1, 2), matrix(1, 2))
  expect_identical(check_net_id_and_initialize(c("id1", "id2"), networks_list), c("id1", "id2"))
  expect_identical(check_net_id_and_initialize(c("id1", "id2", "id3"), list(matrix(1, 2), matrix(1, 2), matrix(1, 2))), c("id1", "id2", "id3"))
})

library("future.apply")
library("future.callr")
library("progressr")
handlers(global = TRUE)
handlers("cli")
# plan(list(tweak("callr", workers = floor(parallelly::availableCores(omit = 2L) / 3L)), tweak("callr", workers = 2L)))
plan("sequential")
options(future.globals.maxSize = Inf)

devtools::load_all()


# data(dorebipartite)
# netlist <- dorebipartite[1:7]
# colsbm_model <- "iid"
# net_id <- NULL
# distribution <- "bernoulli"
# nb_run <- 3L
# global_opts <- list(backend = "no_mc", verbosity = 1L)
# fit_opts <- list()
# fusions_per_step <- 5L
# full_inference <- FALSE

# plan(tweak("callr", workers = parallelly::availableCores(omit = 2L)))
# bbb <- clusterize_bipartite_networks_graphon(
#     netlist,
#     colsbm_model,
#     net_id,
#     distribution,
#     nb_run,
#     global_opts,
#     fit_opts,
#     full_inference = full_inference
# )

nr <- 120L
nc <- 120L

current_model <- "iid"
pi <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
rho <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
eps <- 0.4

alpha_assortative <- matrix(0.3, nrow = 3, ncol = 3) +
  matrix(
    c(
      eps, -0.5 * eps, -0.5 * eps,
      -0.5 * eps, eps, -0.5 * eps,
      -0.5 * eps, -0.5 * eps, eps
    ),
    nrow = 3, byrow = TRUE
  )

alpha_core_periphery <- matrix(0.3, nrow = 3, ncol = 3) +
  matrix(
    c(
      1.5 * eps, eps, 0.5 * eps,
      eps, 0.5 * eps, 0,
      0.5 * eps, 0, -0.5 * eps
    ),
    nrow = 3, byrow = TRUE
  )

alpha_disassortative <- matrix(0.3, nrow = 3, ncol = 3) +
  matrix(
    c(
      -0.5 * eps, eps, 1.5 * eps,
      eps, -0.5 * eps, eps,
      1.5 * eps, eps, -0.5 * eps
    ),
    nrow = 3, byrow = TRUE
  )

assortative_collection <- generate_bipartite_collection(
  nr, nc,
  pi, rho,
  alpha_assortative, 3,
  model = current_model,
  return_memberships = TRUE
)

assortative_incidence <- lapply(
  seq_along(assortative_collection),
  function(m) {
    return(assortative_collection[[m]]$incidence_matrix)
  }
)

assortative_row_clustering <- lapply(
  seq_along(assortative_collection),
  function(m) {
    return(assortative_collection[[m]]$row_clustering)
  }
)

assortative_col_clustering <- lapply(
  seq_along(assortative_collection),
  function(m) {
    return(assortative_collection[[m]]$row_clustering)
  }
)

core_periphery_collection <- generate_bipartite_collection(
  nr, nc,
  pi, rho,
  alpha_core_periphery, 3,
  model = current_model,
  return_memberships = TRUE
)

core_periphery_incidence <- lapply(
  seq_along(core_periphery_collection),
  function(m) {
    return(core_periphery_collection[[m]]$incidence_matrix)
  }
)

core_periphery_row_clustering <- lapply(
  seq_along(core_periphery_collection),
  function(m) {
    return(core_periphery_collection[[m]]$row_clustering)
  }
)

core_periphery_col_clustering <- lapply(
  seq_along(core_periphery_collection),
  function(m) {
    return(core_periphery_collection[[m]]$row_clustering)
  }
)

disassortative_collection <- generate_bipartite_collection(
  nr, nc,
  pi, rho,
  alpha_disassortative, 3,
  model = current_model,
  return_memberships = TRUE
)

disassortative_incidence <- lapply(
  seq_along(disassortative_collection),
  function(m) {
    return(disassortative_collection[[m]]$incidence_matrix)
  }
)

disassortative_row_clustering <- lapply(
  seq_along(disassortative_collection),
  function(m) {
    return(disassortative_collection[[m]]$row_clustering)
  }
)

disassortative_col_clustering <- lapply(
  seq_along(disassortative_collection),
  function(m) {
    return(disassortative_collection[[m]]$row_clustering)
  }
)

real_row_clustering <- append(
  append(
    assortative_row_clustering,
    core_periphery_row_clustering
  ),
  disassortative_row_clustering
)

real_col_clustering <- append(
  append(
    assortative_col_clustering,
    core_periphery_col_clustering
  ),
  disassortative_col_clustering
)

incidence_matrices <- append(
  append(
    assortative_incidence,
    core_periphery_incidence
  ),
  disassortative_incidence
)

netids <- paste0(rep(c("as", "cp", "dis"), each = 3), ".", seq(1, 3))


netlist <- incidence_matrices
net_id <- netids

prefit <- readRDS(here::here("dev", "prefit_9_sim_sep.Rds"))

#
aaa <- clusterize_bipartite_networks_graphon(
  netlist,
  colsbm_model,
  net_id,
  distribution,
  nb_run,
  global_opts,
  fit_opts,
  fit_init = prefit
)

# fits_lbm <- lapply(seq(7, 9), function(id) {
#     sbm::estimateBipartiteSBM(netMat = netlist[[id]], model = "bernoulli", estimOptions = list(verbosity = 1L))
# })

# fit_as <- sbm::estimateBipartiteSBM(netMat = netlist[[1]], model = "bernoulli", estimOptions = list(verbosity = 1L))

# fit_cp <- sbm::estimateBipartiteSBM(netMat = netlist[[4]], model = "bernoulli", estimOptions = list(verbosity = 1L))

# fits_lbm[[4]] <- fit_as
# fits_lbm[[5]] <- fit_cp

# pis <- lapply(fits_lbm, function(fit) {
#     return(fit$blockProp$row)
# })
# rhos <- lapply(fits_lbm, function(fit) {
#     return(fit$blockProp$col)
# })
# alphas <- lapply(fits_lbm, function(fit) {
#     return(fit$connectParam$mean)
# })

# dist_graphon_bipartite_all_permutations(pis = pis[c(1, 2)], rhos = rhos[c(1, 2)], alphas = alphas[c(1, 2)])
# dist_graphon_bipartite_marginals(pis = pis[c(1, 2)], rhos = rhos[c(1, 2)], alphas = alphas[c(1, 2)])

# # 2 and 3
# dist_graphon_bipartite_all_permutations(pis = pis[c(2, 3)], rhos = rhos[c(2, 3)], alphas = alphas[c(2, 3)])
# dist_graphon_bipartite_marginals(pis = pis[c(2, 3)], rhos = rhos[c(2, 3)], alphas = alphas[c(2, 3)])

# # 1 and 3
# dist_graphon_bipartite_all_permutations(pis = pis[c(1, 3)], rhos = rhos[c(1, 3)], alphas = alphas[c(1, 3)])
# dist_graphon_bipartite_marginals(pis = pis[c(1, 3)], rhos = rhos[c(1, 3)], alphas = alphas[c(1, 3)])

# # 1 and as
# dist_graphon_bipartite_all_permutations(pis = pis[c(1, 4)], rhos = rhos[c(1, 4)], alphas = alphas[c(1, 4)])
# dist_graphon_bipartite_marginals(pis = pis[c(1, 4)], rhos = rhos[c(1, 4)], alphas = alphas[c(1, 4)])

# # 2 and as
# dist_graphon_bipartite_all_permutations(pis = pis[c(2, 4)], rhos = rhos[c(2, 4)], alphas = alphas[c(2, 4)])
# dist_graphon_bipartite_marginals(pis = pis[c(2, 4)], rhos = rhos[c(2, 4)], alphas = alphas[c(2, 4)])

# # 3 and as
# dist_graphon_bipartite_all_permutations(pis = pis[c(3, 4)], rhos = rhos[c(3, 4)], alphas = alphas[c(3, 4)])
# dist_graphon_bipartite_marginals(pis = pis[c(3, 4)], rhos = rhos[c(3, 4)], alphas = alphas[c(3, 4)])

# # 1 and cp
# dist_graphon_bipartite_all_permutations(pis = pis[c(1, 5)], rhos = rhos[c(1, 5)], alphas = alphas[c(1, 5)])
# dist_graphon_bipartite_marginals(pis = pis[c(1, 5)], rhos = rhos[c(1, 5)], alphas = alphas[c(1, 5)])

# # 2 and cp
# dist_graphon_bipartite_all_permutations(pis = pis[c(2, 5)], rhos = rhos[c(2, 5)], alphas = alphas[c(2, 5)])
# dist_graphon_bipartite_marginals(pis = pis[c(2, 5)], rhos = rhos[c(2, 5)], alphas = alphas[c(2, 5)])

# # 3 and cp
# dist_graphon_bipartite_all_permutations(pis = pis[c(3, 5)], rhos = rhos[c(3, 5)], alphas = alphas[c(3, 5)])
# dist_graphon_bipartite_marginals(pis = pis[c(3, 5)], rhos = rhos[c(3, 5)], alphas = alphas[c(3, 5)])
# aaa <- clusterize_bipartite_networks_graphon(
#     incidence_matrices,
#     colsbm_model = "iid",
#     net_id = netids,
#     distribution,
#     nb_run,
#     global_opts,
#     fit_opts
# )

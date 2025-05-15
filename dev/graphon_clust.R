library("future.apply")
library("future.callr")
library("progressr")
handlers(global = TRUE)
handlers("cli")
# plan(list(tweak("callr", workers = floor(parallelly::availableCores(omit = 2L) / 3L)), tweak("callr", workers = 2L)))
plan("sequential")
options(future.globals.maxSize = Inf)

devtools::load_all()


data(dorebipartite)
netlist <- dorebipartite[1:7]
colsbm_model <- "iid"
net_id <- NULL
distribution <- "bernoulli"
nb_run <- 1L
global_opts <- list(backend = "no_mc", verbosity = 1L)
fit_opts <- list()
fusions_per_step <- 5L
full_inference <- FALSE

plan(tweak("callr", workers = parallelly::availableCores(omit = 2L)))
# classic_penalty <- clusterize_bipartite_networks(
#     netlist,
#     colsbm_model,
#     net_id,
#     distribution,
#     nb_run,
#     global_opts,
#     fit_opts,
#     full_inference = full_inference
# )

no_penalty <- clusterize_bipartite_networks(
    netlist,
    colsbm_model,
    net_id,
    distribution,
    nb_run,
    global_opts,
    fit_opts = list(penalty_factor = 0),
    full_inference = full_inference
)

small_penalty <- clusterize_bipartite_networks(
    netlist,
    colsbm_model,
    net_id,
    distribution,
    nb_run,
    global_opts,
    fit_opts = list(penalty_factor = 0.1),
    full_inference = full_inference
)

penalty_005 <- clusterize_bipartite_networks(
    netlist,
    colsbm_model,
    net_id,
    distribution,
    nb_run,
    global_opts,
    fit_opts = list(penalty_factor = 0.05),
    full_inference = full_inference
)

large_penalty <- clusterize_bipartite_networks(
    netlist,
    colsbm_model,
    net_id,
    distribution,
    nb_run,
    global_opts,
    fit_opts = list(penalty_factor = 10),
    full_inference = full_inference
)

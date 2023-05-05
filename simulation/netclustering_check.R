require("ggplot2")
require("tictoc")

devtools::load_all("R/")

# Generation of data
incidence_matrices <- NULL

if (!exists("arg")) {
    arg <- commandArgs(trailingOnly = TRUE)
}

if (identical(arg, character(0))) {
    number_of_net <- length(incidence_matrices)
} else {
    number_of_net <- as.numeric(arg)
}

tic()
list_collection <- clusterize_bipartite_networks(
    netlist = incidence_matrices[1:number_of_net],
    colsbm_model = "pirho",
    global_opts = list(
        nb_cores = parallel::detectCores() - 1, verbosity = 1,
        plot_details = 0
    ),
    full_inference = TRUE
)
toc()

saveRDS(list_collection, file = paste0(
    "simulation/data/", 
    "simulated_collection_clustering",
    number_of_net, "-",
    format(Sys.time(), "%d-%m-%y-%X"),
    ".Rds"
))

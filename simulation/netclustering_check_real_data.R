require("ggplot2")
require("tictoc")
require("tidyverse")

devtools::load_all("R/")

# Importation of data
if (!file.exists("simulation/data/dore-matrices.Rds")) {
    interaction_data <- read.table(file = "simulation/data/interaction-data.txt", sep = "\t", header = TRUE)

    seq_ids_network_aggreg <- unique(interaction_data$id_network_aggreg)
    names_aggreg_networks <- sapply(
        seq_ids_network_aggreg,
        function(id) {
            paste0(
                unique(interaction_data[which(interaction_data$id_network_aggreg == id), ]$web),
                collapse = "+"
            )
        }
    )
    # Computation of incidence matrices
    incidence_matrices <- lapply(
        seq_ids_network_aggreg,
        function(m) {
            current_interaction_data <- interaction_data[which(interaction_data$id_network_aggreg == m), ] %>%
                mutate(
                    plantaggreg = paste(plantorder,
                        plantfamily, plantgenus, plantspecies,
                        sep = "-"
                    ),
                    insectaggreg = paste(insectorder,
                        insectfamily, insectgenus, insectspecies,
                        sep = "-"
                    )
                )
            current_interaction_data <- table(current_interaction_data$plantaggreg, current_interaction_data$insectaggreg)

            current_incidence_matrix <- matrix(current_interaction_data,
                ncol = ncol(current_interaction_data), dimnames = dimnames(current_interaction_data)
            )

            current_incidence_matrix[which(current_incidence_matrix > 0)] <- 1
            return(current_incidence_matrix)
        }
    )

    names(incidence_matrices) <- names_aggreg_networks

    saveRDS(incidence_matrices, file = "simulation/data/dore-matrices.Rds")
} else {
    incidence_matrices <- readRDS(file = "simulation/data/dore-matrices.Rds")
}

if (!exists("arg")) {
    arg <- commandArgs(trailingOnly = TRUE)
}

if (identical(arg, character(0))) {
    number_of_net <- length(incidence_matrices)
    model <- "pirho"
    nb_run <- 3
} else {
    number_of_net <- as.numeric(arg[1])
    model <- arg[2]
    nb_run <- as.numeric(arg[3])
}
print(number_of_net)
print(model)
print(nb_run)
tic()
list_collection <- clusterize_bipartite_networks(
    netlist = incidence_matrices[1:number_of_net],
    net_id = names(incidence_matrices)[1:number_of_net],
    colsbm_model = model,
    nb_run = nb_run,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1, verbosity = 4,
        plot_details = 0
    ),
    silent_parallelization = TRUE
)
toc()

saveRDS(list_collection, file = paste0(
    "simulation/data/",
    "dore_collection_clustering_nb_run", nb_run,"_",model,"_",
    number_of_net, "networks_",
    format(Sys.time(), "%d-%m-%y-%X"),
    ".Rds"
))

require("colSBM")
require("ggplot2")
require("ggnewscale")

devtools::load_all("R/")

# Importation of data
if (!file.exists("simulation/data/dore-matrices.Rds")) {
    interaction_data <- read.table(file = "simulation/data/interaction-data.txt", sep = "\t", header = TRUE)

    seq_ids_network_aggreg <- unique(interaction_data$id_network_aggreg)

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

    names(incidence_matrices) <- seq_ids_network_aggreg

    saveRDS(incidence_matrices, file = "simulation/data/dore-matrices.Rds")
} else {
    incidence_matrices <- readRDS(file = "simulation/data/dore-matrices.Rds")
}

list_collection <- clusterize_bipartite_networks(
    netlist = incidence_matrices[1:3],
    colsbm_model = "pirho",
    global_opts = list(
        nb_cores = parallel::detectCores() - 1, verbosity = 0,
        plot_details = 0
    )
)

saveRDS(list_collection, file = paste0(
    "simulations/data/", 
    "dore-collection-clustering", 
    Sys.time(),
    ".Rds"
))
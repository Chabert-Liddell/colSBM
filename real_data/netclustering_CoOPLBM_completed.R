require("ggplot2")
require("tictoc")
require("tidyverse")

devtools::load_all("R/")

# Importation of data
interaction_data <- read.table(file = "real_data/data/interaction-data.txt", sep = "\t", header = TRUE)

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

# Full data
if (!file.exists("real_data/data/dore-matrices.Rds")) {
  # Computation of incidence matrices
  incidence_matrices <- lapply(
    seq_ids_network_aggreg,
    function(m) {
      current_interaction_data <- interaction_data[which(interaction_data$id_network_aggreg == m), ] %>%
        mutate(
          plantaggreg = paste(plantorder,
            plantfamily, plantgenus, plantspecies,
            sep = " "
          ),
          insectaggreg = paste(insectorder,
            insectfamily, insectgenus, insectspecies,
            sep = " "
          )
        )
      current_interaction_data <- table(current_interaction_data$insectaggreg, current_interaction_data$plantaggreg)

      current_incidence_matrix <- matrix(current_interaction_data,
        ncol = ncol(current_interaction_data), dimnames = dimnames(current_interaction_data)
      )

      current_incidence_matrix[which(current_incidence_matrix > 0)] <- 1
      return(current_incidence_matrix)
    }
  )

  names(incidence_matrices) <- names_aggreg_networks

  saveRDS(incidence_matrices, file = "real_data/data/dore-matrices.Rds")
} else {
  incidence_matrices <- readRDS(file = "real_data/data/dore-matrices.Rds")
}
# Emre completed Data
completed_networks_ids <- as.numeric(names(readRDS("real_data/data/Data.rds")))
completed_networks_names <- sapply(
  completed_networks_ids,
  function(id) {
    paste0(
      unique(interaction_data[which(interaction_data$id_network_aggreg == id), ]$web),
      collapse = "+"
    )
  }
)

uncompleted <- incidence_matrices[match(completed_networks_names, names(incidence_matrices))]
point_2_completed <- setNames(readRDS("real_data/data/completed0.2.rds"), completed_networks_names)
point_5_completed <- setNames(readRDS("real_data/data/completed0.5.rds"), completed_networks_names)
random_completed <- setNames(readRDS("real_data/data/completedrandom.rds"), completed_networks_names)

if (!exists("arg")) {
  arg <- commandArgs(trailingOnly = TRUE)
}

if (identical(arg, character(0))) {
  number_of_net <- length(point_2_completed)
  model <- "pirho"
  nb_run <- 3
  source_data <- eval(as.name("point_2_completed"))
  data_name <- "point_2_completed"
} else {
  number_of_net <- as.numeric(arg[1])
  model <- arg[2]
  nb_run <- as.numeric(arg[3])
  source_data <- eval(as.name(arg[4]))
  data_name <- arg[4]
}
print(number_of_net)
print(model)
print(nb_run)
print(data_name)
tic()
list_collection <- clusterize_bipartite_networks(
  netlist = source_data[1:number_of_net],
  net_id = names(source_data)[1:number_of_net],
  colsbm_model = model,
  nb_run = nb_run,
  global_opts = list(
    nb_cores = parallel::detectCores() - 1, verbosity = 4,
    plot_details = 0
  ),
  silent_parallelization = FALSE
)
toc()

saveRDS(list_collection, file = paste0(
  "real_data/data/",
  "dore_",
  data_name,
  "_collection_clustering_nb_run_", nb_run, "_", model, "_",
  number_of_net, "_networks_",
  format(Sys.time(), "%d-%m-%y-%X"),
  ".Rds"
))

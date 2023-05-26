require("ggplot2")
require("tictoc")

devtools::load_all("R/")

iid_clustering <- readRDS("simulation/data/dore_collection_clustering_nb_run1_iid_123networks_24-05-23-21:40:42.Rds")

iid_best_partition <- extract_bipartite_best_partition(iid_clustering)

rho_clustering <- readRDS("simulation/data/dore_collection_clustering_nb_run1_rho_123networks_25-05-23-13:58:30.Rds")
rho_best_partition <- extract_bipartite_best_partition(rho_clustering)

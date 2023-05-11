require("ggplot2")
require("tictoc")

devtools::load_all("R/")

result_clustering <- readRDS("simulation/data/simulated_collection_clustering_rho_10-05-23-14:40:46.Rds")

list_clustering <- lapply(
    seq_along(result_clustering), function(s) result_clustering[[s]]$list_of_clusterings
)

list_best_partition <- lapply(
    seq_along(list_clustering), function(s) {
        list(
            epsilon = result_clustering[[s]]$epsilon,
            best_partition = extract_bipartite_best_partition(list_clustering[[s]])
        )
    }
)

require("ggplot2")
devtools::load_all("R/")

## Preparing data for analysis
interaction_data <- read.table(file = "real_data/data/interaction-data.txt", sep = "\t", header = TRUE)

insect_orders <- unique(interaction_data$insectorder)
plant_family <- unique(interaction_data$plantfamily)

# All results
iid_clustering <- readRDS("real_data/data/dore_collection_clustering_nb_run1_iid_123networks_24-05-23-21:40:42.Rds")
iid_best_partition <- extract_best_bipartite_partition(iid_clustering)
iid_unlist <- unlist(iid_best_partition)


rho_clustering <- readRDS("real_data/data/dore_collection_clustering_nb_run1_rho_123networks_25-05-23-13:58:30.Rds")
rho_best_partition <- extract_best_bipartite_partition(rho_clustering)
rho_unlist <- unlist(rho_best_partition)

pi_clustering <- readRDS("real_data/data/dore_collection_clustering_nb_run1_pi_123networks_25-05-23-17:31:25.Rds")
pi_best_partition <- extract_best_bipartite_partition(pi_clustering)
pi_unlist <- unlist(pi_best_partition)

pirho_clustering <- readRDS("real_data/data/dore_collection_clustering_nb_run1_pirho_123networks_26-05-23-19:22:55.Rds")
pirho_best_partition <- extract_best_bipartite_partition(pirho_clustering)
pirho_unlist <- unlist(pirho_best_partition)

### Matching taxonomy
taxonomy_in_clusters <- function(unlisted_model) {
  lapply(seq_len(length(unlisted_model)), function(col_idx) {
    # Per collection
    # Empty init
    insect_count <- t(sapply(insect_orders, function(order) {
      out_count <- rep(0, unlisted_model[[col_idx]]$Q[2])
      out_count
    }))

    plant_count <- t(sapply(plant_family, function(order) {
      out_count <- rep(0, unlisted_model[[col_idx]]$Q[1])
      out_count
    }))

    for (m in seq.int(unlisted_model[[col_idx]]$M)) {
      #### Insect
      insect_names <- names(unlisted_model[[col_idx]]$Z[[1]][[2]])

      insect_count <- insect_count + t(sapply(insect_orders, function(order) {
        out_count <- rep(0, unlisted_model[[col_idx]]$Q[2])
        names(out_count) <- seq.int(unlisted_model[[col_idx]]$Q[2])
        insect_count <- table(unlisted_model[[col_idx]]$Z[[m]][[2]][grep(order, insect_names)])
        out_count[names(insect_count)] <- insect_count
        out_count
      }))
      #### Plants
      plant_names <- names(unlisted_model[[col_idx]]$Z[[1]][[1]])

      plant_count <- t(sapply(plant_family, function(order) {
        out_count <- rep(0, unlisted_model[[col_idx]]$Q[1])
        names(out_count) <- seq.int(unlisted_model[[col_idx]]$Q[1])
        plant_count <- table(unlisted_model[[col_idx]]$Z[[m]][[1]][grep(order, plant_names)])
        out_count[names(plant_count)] <- plant_count
        out_count
      }))
    }
    return(list(insects = insect_count, plants = plant_count))
  })
}

taxonomy_remove_empty <- function(taxonomy_collections_list) {
  lapply(taxonomy_collections_list, function(collection) {
    list(
      insects = collection$insects[which(rowSums(collection$insects != 0) > 0), ],
      plants = collection$plants[which(rowSums(collection$plants != 0) > 0), ]
    )
  })
}

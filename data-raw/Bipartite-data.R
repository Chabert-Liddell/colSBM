require("dplyr")

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

# Computation of incidence matrices
incidence_matrices <- lapply(
  seq_ids_network_aggreg,
  function(m) {
    current_interaction_data <- interaction_data[which(interaction_data$id_network_aggreg == m), ] |>
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

dorebipartite <- incidence_matrices[1:15]

usethis::use_data(dorebipartite)

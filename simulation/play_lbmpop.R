# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")
source("R/R6class-lbmpop.R")

set.seed(1234)
eps <- 0.05
M <- 3

# Cas hote parasite
# Modulaire : avec les diagonaux qui sont x3 et x4 des petites probas non diag

nr_modular <- 200
nr_modular <- c(1 / 10 * nr_modular, 1 / 5 * nr_modular, nr_modular)
nc_modular <- 500
nc_modular <- c(1 / 10 * nc_modular, 1 / 5 * nc_modular, nc_modular)

pir_modular <- c(0.6,0.4)
pic_modular <- c(0.6,0.4)

Q_modular <- c(length(pir_modular), length(pic_modular))

alpha_modular <- matrix(
    c(
        15*eps, eps,
        eps, 10*eps
    ),
    nrow = Q_modular[1], ncol = Q_modular[2], byrow = TRUE
)

modular_bipartite_collection <- generate_bipartite_collection(nr_modular, nc_modular, pir_modular, pic_modular, alpha_modular, M)

modular_bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    modular_bipartite_collection[[m]]$incidence_matrix
})

### Cas plante-pollinisateurs double emboitements avec (1,3) partagé entre les emboîtements
nr_nested <- 200
nr_nested <- c(1/10*nr_nested, 1/5*nr_nested, nr_nested)
nc_nested <- 500
nc_nested <- c(1/10*nc_nested, 1/5*nc_nested, nc_nested)

pir_nested <- c(1)
pic_nested <- c(1)

Q_nested <- c(length(pir_nested), length(pic_nested))

alpha_nested <- matrix(
    c(
        0.25
    ),
    nrow = Q_nested[1], ncol = Q_nested[2], byrow = TRUE
)

nested_bipartite_collection <- generate_bipartite_collection(nr_nested, nc_nested, pir_nested, pic_nested, alpha_nested, M)

nested_bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    nested_bipartite_collection[[m]]$incidence_matrix
})

# La fusion des deux précédents

bipartite_collection <- list(modular_bipartite_collection, nested_bipartite_collection)
# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq_along(bipartite_collection), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

# Z <- lapply(seq_along(bipartite_collection), function(m) {
#     list(
#         bipartite_collection[[m]]$row_clustering,
#         bipartite_collection[[m]]$col_clustering
#     )
# })


# # mylbmpop$burn_in()
# # mylbmpop$moving_window(mylbmpop$best_fit$Q)
# # mylbmpop$compute_sep_LBM_ICL()
# mylbmpop$optimize()

# colLBM models
cat("\n Hote-Parasite (Modulaire) :\n")

modular_lbmpop <- lbmpop$new(
    netlist = modular_bipartite_collection_incidence,
    free_density = FALSE,
    free_mixture = FALSE,
    global_opts = list(verbosity = 4, plot_details = 1, nb_cores = 2)
)

modular_lbmpop$optimize()

cat("\n Plantes-Pollinisateurs (Double emboîté) :\n")

nested_lbmpop <- lbmpop$new(
    netlist = nested_bipartite_collection_incidence,
    free_density = FALSE,
    free_mixture = FALSE,
    global_opts = list(verbosity = 4, plot_details = 1)
)

nested_lbmpop$optimize()

cat("\n La collection associant les deux précédentes:\n")

bad_lbmpop <- lbmpop$new(
    netlist = bipartite_collection_incidence,
    free_density = FALSE,
    free_mixture = FALSE,
    global_opts = list(
        verbosity = 4,
        plot_details = 1,
        Q1_max = 10,
        Q2_max = 10
    )
)

bad_lbmpop$optimize()

modular <- estimateBipartiteSBM(bipartite_collection_incidence[[1]])
nested <- estimateBipartiteSBM(bipartite_collection_incidence[[2]])

beepr::beep(5)

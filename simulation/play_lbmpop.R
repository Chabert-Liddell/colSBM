# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

# require("visNetwork", quietly = T)
# require("igraph", quietly = T)
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")
source("R/R6class-lbmpop.R")

set.seed(1234)
eps <- 0.05
M <- 3

# Cas hote parasite
# Modulaire : avec les diagonaux qui sont x3 et x4 des petites probas non diag

nr_hp <- 200
nr_hp <- c(1 / 10 * nr_hp, 1 / 5 * nr_hp, nr_hp)
nc_hp <- 500
nc_hp <- c(1 / 10 * nc_hp, 1 / 5 * nc_hp, nc_hp)

pir_hp <- c(0.6,0.4)
pic_hp <- c(0.6,0.1,0.3)

Q_hp <- c(length(pir_hp), length(pic_hp))

alpha_hp <- matrix(
    c(
        15*eps, eps, eps,
        eps, 10*eps, eps,
        eps, eps, 7*eps
    ),
    nrow = Q_hp[1], ncol = Q_hp[2], byrow = TRUE
)

hp_bipartite_collection <- generate_bipartite_collection(nr_hp, nc_hp, pir_hp, pic_hp, alpha_hp, M)

hp_bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    hp_bipartite_collection[[m]]$incidence_matrix
})

### Cas plante-pollinisateurs double emboitements avec (1,3) partagé entre les emboîtements
nr_pp <- 200
nr_pp <- c(1/10*nr_pp, 1/5*nr_pp, nr_pp)
nc_pp <- 500
nc_pp <- c(1/10*nc_pp, 1/5*nc_pp, nc_pp)

pir_pp <- c(0.4, 0.3, 0.2, 0.1)
pic_pp <- c(0.4, 0.3, 0.2, 0.1)

Q_pp <- c(length(pir_pp), length(pic_pp))

alpha_pp <- matrix(
    c(
        15*eps, 10*eps, 10*eps, eps,
        10*eps, eps, eps, eps,
        eps, eps, 15*eps, 10*eps,
        eps, eps, 10*eps, eps
    ),
    nrow = Q_pp[1], ncol = Q_pp[2], byrow = TRUE
)

pp_bipartite_collection <- generate_bipartite_collection(nr_pp, nc_pp, pir_pp, pic_pp, alpha_pp, M)

pp_bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    pp_bipartite_collection[[m]]$incidence_matrix
})

# La fusion des deux précédents

bipartite_collection <- append(hp_bipartite_collection, pp_bipartite_collection)
# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq_along(bipartite_collection), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq_along(bipartite_collection), function(m) {
    list(
        bipartite_collection[[m]]$row_clustering,
        bipartite_collection[[m]]$col_clustering
    )
})


# # mylbmpop$burn_in()
# # mylbmpop$moving_window(mylbmpop$best_fit$Q)
# # mylbmpop$compute_sep_LBM_ICL()
# mylbmpop$optimize()

# colLBM models
cat("\n Hote-Parasite (Modulaire) :\n")

hp_lbmpop <- lbmpop$new(
    netlist = hp_bipartite_collection_incidence,
    free_density = FALSE,
    free_mixture = FALSE,
    global_opts = list(verbosity = 4, plot_details = 1)
)

hp_lbmpop$optimize()

cat("\n Plantes-Pollinisateurs (Double emboîté) :\n")

pp_lbmpop <- lbmpop$new(
    netlist = pp_bipartite_collection_incidence,
    free_density = FALSE,
    free_mixture = FALSE,
    global_opts = list(verbosity = 4, plot_details = 1)
)

pp_lbmpop$optimize()

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
beepr::beep(3)
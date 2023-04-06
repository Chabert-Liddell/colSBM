# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

devtools::load_all(path = "R/")

set.seed(1234)

verbose <- TRUE
test_alea <- TRUE

eps <- 0.05
M <- 3
nr <- 100
nc <- 250

pir <- c(0.2, 0.8)
pic <- c(0.2, 0.8)

Q <- c(length(pir), length(pic))

alpha <- matrix(
    c(
        0.9, eps,
        eps, 0.8
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})

choosed_bisbmpop <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence, 
    colsbm_model = "iid", 
    global_opts = list(nb_cores = 3)
)

ari_sums <- sapply(
    seq_along(choosed_bisbmpop$best_fit$Z),
    function(m) {
        sum(c(
            aricode::ARI(
                Z[[m]][[1]],
                choosed_bisbmpop$best_fit$Z[[m]][[1]]
            ),
            aricode::ARI(
                Z[[m]][[2]],
                choosed_bisbmpop$best_fit$Z[[m]][[2]]
            )
        ))
    }
)
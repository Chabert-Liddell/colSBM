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
M <- 2
nr <- 100
nc <- 250

pir1 <- c(0.2, 0, 0.8)
pir2 <- c(0.4, 0.6, 0)

pic <- c(0.2, 0.8)


Q <- c(length(pir1), length(pic))

# Make a non common alpha structure

alpha <- matrix(
    c(
        0.4, eps,
        eps, 0.5,
        0.2, eps
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- list(
    generate_bipartite_network(nr, nc, pir1, pic, alpha),
    generate_bipartite_network(nr, nc, pir2, pic, alpha)#,
    #generate_bipartite_network(nr, nc, pir3, pic, alpha)
)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})


mybisbmpop <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence, 
    colsbm_model = "pi", 
    global_opts = list(nb_cores = 3)
)

# ari_sums <- sapply(
#     seq_along(choosed_bisbmpop$best_fit$Z),
#     function(m) {
#         sum(c(
#             aricode::ARI(
#                 Z[[m]][[1]],
#                 choosed_bisbmpop$best_fit$Z[[m]][[1]]
#             ),
#             aricode::ARI(
#                 Z[[m]][[2]],
#                 choosed_bisbmpop$best_fit$Z[[m]][[2]]
#             )
#         ))
#     }
# )
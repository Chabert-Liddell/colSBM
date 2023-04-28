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

pic1 <- c(0.2, 0, 0.8)
pic2 <- c(0.4, 0.6, 0)

pir <- c(0.2, 0.8)


Q <- c(length(pir), length(pic1))

# Make a non common alpha structure

alpha <- matrix(
    c(
        0.4, eps, eps,
        eps, 0.5, eps
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- list(
    generate_bipartite_network(nr, nc, pir, pic1, alpha),
    generate_bipartite_network(nr, nc, pir, pic2, alpha)
)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})

mybisbmpop <- bisbmpop$new(
    netlist = bipartite_collection_incidence,
    free_mixture_col = TRUE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 4,
        parallelization_vector = c(FALSE, TRUE)
    )
)

mybisbmpop$optimize()

# choosed_bisbmpop <- estimate_colBiSBM(
#     netlist = bipartite_collection_incidence, 
#     colsbm_model = "iid", 
#     global_opts = list(nb_cores = 3)
# )

ari_sums <- sapply(
    seq_along(mybisbmpop$best_fit$Z),
    function(m) {
        c(
            aricode::ARI(
                Z[[m]][[1]],
                mybisbmpop$best_fit$Z[[m]][[1]]
            ),
            aricode::ARI(
                Z[[m]][[2]],
                mybisbmpop$best_fit$Z[[m]][[2]]
            )
        )
    }
)

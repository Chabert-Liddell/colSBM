# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

devtools::load_all(path = "R/")

tic()

eps <- 0.05
nr <- 250
nc <- 250

pir1 <- c(0.4, 0, 0.6)
pir2 <- c(0.4, 0.6, 0)
pir3 <- c(0, 0.4, 0.6)

pic1 <- c(0.6, 0, 0.4)
pic2 <- c(0.6, 0.4, 0)
pic3 <- c(0, 0.6, 0.4)

Q <- c(length(pir1), length(pic1))

# Make a non common alpha structure

alpha <- matrix(
    c( # 12   2    1
        0.6, 0.25, eps,# 12
        eps, 0.8, 0,# 2
        0.2, 0, 0.4 # 1NB
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- list(
    generate_bipartite_network(nr, nc, pir1, pic1, alpha),
    generate_bipartite_network(nr, nc, pir2, pic2, alpha),
    generate_bipartite_network(nr, nc, pir3, pic3, alpha),
    generate_bipartite_network(nr, nc, pir2, pic3, alpha),
    generate_bipartite_network(nr, nc, pir3, pic2, alpha),
    generate_bipartite_network(nr, nc, pir1, pic2, alpha),
    generate_bipartite_network(nr, nc, pir2, pic1, alpha)
)

M <- length(bipartite_collection)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
})


## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})

row_clusterings <- lapply(seq_along(bipartite_collection), function(m) {
    return(bipartite_collection[[m]]$row_clustering)
})

col_clusterings <- lapply(seq_along(bipartite_collection), function(m) {
    return(bipartite_collection[[m]]$col_clustering)
})

full_row_clustering <- as.vector(sapply(
    seq.int(M),
    function(m) row_clusterings[[m]]
))

full_col_clustering <- as.vector(sapply(
    seq.int(M),
    function(m) col_clusterings[[m]]
))

pi <- list(
    list(pir1, pic1),
    list(pir2, pic2),
    list(pir3, pic3),
    list(pir2, pic3),
    list(pir3, pic2),
    list(pir1, pic2),
    list(pir2, pic1)
)

Cpi <- vector(mode = "list", length = 2)
Cpi[[1]] <- vapply(seq(M), function(m) {
    pi[[m]][[1]] > 0
},
FUN.VALUE = rep(TRUE, Q[1])
)
Cpi[[2]] <- vapply(seq(M), function(m) {
    pi[[m]][[2]] > 0
},
FUN.VALUE = rep(TRUE, Q[2])
)

Calpha <- tcrossprod(Cpi[[1]], Cpi[[2]]) > 0

mybisbmpop <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "iid",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 1,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop1 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "pi",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1, 
        verbosity = 1,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop2 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "rho",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 1,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop3 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "pirho",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 1,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

compute_double_ARI <- function(model) {
    model_row_Z <- as.vector(sapply(
        seq.int(model$M),
        function(m) model$Z[[m]][[1]]
    ))

    model_col_Z <- as.vector(sapply(
        seq.int(model$M),
        function(m) model$Z[[m]][[2]]
    ))

    return(list(
        aricode::ARI(model_row_Z, full_row_clustering),
        aricode::ARI(model_col_Z, full_col_clustering)
    ))
}

ari_sums <- sapply(
    seq_along(mybisbmpop3$best_fit$Z),
    function(m) {
        c(
            aricode::ARI(
                Z[[m]][[1]],
                mybisbmpop3$best_fit$Z[[m]][[1]]
            ),
            aricode::ARI(
                Z[[m]][[2]],
                mybisbmpop3$best_fit$Z[[m]][[2]]
            )
        )
    }
)
toc()
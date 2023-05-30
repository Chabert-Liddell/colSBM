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

pir1 <- c(0.4, 0, 0.5, 0.1)
pir2 <- c(0.4, 0.5, 0, 0.1)

pic1 <- c(0.4, 0, 0.5, 0.1)
pic2 <- c(0.4, 0.5, 0, 0.1)


Q <- c(length(pir1), length(pic1))

# Make a non common alpha structure

alpha <- matrix(
    c( # 12   2    1
        0.6, 0.25, eps, 0.7, # 12
        eps, 0.8, 0, 0.2,# 2
        0.2, 0, 0.4, 0.45,# 1NB
        eps, 0.3, 0.1, 0.7
    ), nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- list(
    generate_bipartite_network(nr, nc, pir1, pic1, alpha),
    generate_bipartite_network(nr, nc, pir2, pic2, alpha)
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
    list(pir2, pic2)
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

# mybisbmpop <- estimate_colBiSBM(
#     netlist = bipartite_collection_incidence,
#     colsbm_model = "iid",
#     silent_parallelization = FALSE,
#     global_opts = list(
#         nb_cores = parallel::detectCores() - 1,
#         verbosity = 4,
#         parallelization_vector = c(TRUE, TRUE)
#     )
# )

# mybisbmpop1 <- estimate_colBiSBM(
#     netlist = bipartite_collection_incidence,
#     colsbm_model = "pi",
#     silent_parallelization = FALSE,
#     global_opts = list(
#         nb_cores = parallel::detectCores() - 1, 
#         verbosity = 4,
#         parallelization_vector = c(TRUE, TRUE)
#     )
# )

# mybisbmpop2 <- estimate_colBiSBM(
#     netlist = bipartite_collection_incidence,
#     colsbm_model = "rho",
#     silent_parallelization = FALSE,
#     global_opts = list(
#         nb_cores = parallel::detectCores() - 1,
#         verbosity = 4,
#         parallelization_vector = c(TRUE, TRUE)
#     )
# )

mybisbmpop3 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "pirho",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 4,
        parallelization_vector = c(TRUE, TRUE)
    )
)
toc()
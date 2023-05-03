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

pir1 <- c(0.2, 0, 0.8)
pir2 <- c(0.4, 0.6, 0)

pic1 <- c(0.6, 0, 0.4)
pic2 <- c(0.4, 0.6, 0)

Q <- c(length(pir1), length(pic1))

# Make a non common alpha structure

alpha <- matrix(
    c( # 12   2    1
        0.6, eps, eps,# 12
        eps, 0.8, 0,# 2
        0.2, 0, 0.4 # 1NB
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

pi <- list(list(pir1, pic1), list(pir2, pic2))

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
        verbosity = 4,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop1 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "pi",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1, 
        verbosity = 4,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop2 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "rho",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 4,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

mybisbmpop3 <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence,
    colsbm_model = "pirho",
    silent_parallelization = FALSE,
    global_opts = list(
        nb_cores = parallel::detectCores() - 1,
        verbosity = 4,
        parallelization_vector = c(TRUE, TRUE, FALSE)
    )
)

pirho_best_fit <- mybisbmpop3$best_fit$clone()
pirho_best_fit$MAP$Z <- Z
pirho_best_fit$MAP$pi <- pi
pirho_best_fit$Cpi <- Cpi
pirho_best_fit$Calpha <- Calpha
pirho_best_fit$compute_vbound(MAP = TRUE)
pirho_best_fit$compute_penalty(MAP = TRUE)
pirho_best_fit$compute_BICL(MAP = TRUE)


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
toc()
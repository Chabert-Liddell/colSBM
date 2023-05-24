# Sourcing all necessary files
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

devtools::load_all(path = "R/")

set.seed(1234)

eps <- 0.05
nr <- 100
nc <- 250

pir1 <- c(0.2, 0.8)
pir2 <- c(0.4, 0.6)
pir3 <- c(0.3, 0.7)

pic <- c(0.2, 0.8)


Q <- c(length(pir1), length(pic))


alpha <- matrix(
    c(
        0.9, eps,
        eps, 0.8
    ),
    nrow = Q[1], ncol = Q[2], byrow = TRUE
)
results <- bettermc::mclapply(seq.int(5), function(r) {
    bipartite_collection <- rep(list(
        generate_bipartite_network(nr, nc, pir1, pic, alpha),
        generate_bipartite_network(nr, nc, pir2, pic, alpha),
        generate_bipartite_network(nr, nc, pir3, pic, alpha)
    ), r)
    M <- length(bipartite_collection)
    # This is a list of the M incidence matrices
    bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
        bipartite_collection[[m]]$incidence_matrix
    })
    start_time <- Sys.time()
    estimate_colBiSBM(
        netlist = bipartite_collection_incidence,
        colsbm_model = "pi",
        global_opts = list(
            nb_cores = parallel::detectCores() - 1,
            verbosity = 0
        ),
        silent_parallelization = TRUE
    )
    end_time <- Sys.time()
    cat("\nFinished", r)
    return(c(M, end_time - start_time))
}, mc.cores = parallel::detectCores() - 1, mc.progress = TRUE)

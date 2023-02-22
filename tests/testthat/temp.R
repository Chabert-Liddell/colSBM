# Sourcing all necessary files
require("sbm", quietly = T)
require("aricode", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
# require("visNetwork", quietly = T)
# require("igraph", quietly = T)
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

eps <- 0.05
M <- 5
nr <- 100
nc <- 250

pir <- c(0.2, 0.8)
pic <- c(0.2, 0.3, 0.5)

Q <- c(length(pir), length(pic))

alpha <- matrix(
    c(
        0.9, eps, eps,
        eps, 0.8, eps
    ),
    nrow = Q[1], ncol = Q[2], byrow = TRUE
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


fitColHCA.minibatch <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    init_method = "hca",
    fit_opts = list(verbosity = 0, minibatch = TRUE)
)

fitColHCA.minibatch$optimize()

fitColHCA.vanilla <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    init_method = "hca",
    fit_opts = list(verbosity = 0, minibatch = FALSE)
)

fitColHCA.vanilla$optimize()
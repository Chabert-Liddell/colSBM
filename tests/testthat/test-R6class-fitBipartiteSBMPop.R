# Sourcing all necessary files
require("sbm")
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

set.seed(1234)

# Generate a bipartite collection fixed nr and nc

eps <- 0.05

nr <- 10
nc <- 25

pir <- c(0.1, 0.2, 0.4, 0.3)
pic <- c(0.7, 0.1, 0.05, 0.05, 0.1)

alpha <- matrix(
    c(
        0.9, eps, eps, eps, eps,
        eps, 0.9, eps, eps, eps,
        eps, eps, 0.9, eps, eps,
        eps, eps, eps, 0.9, eps
    ),
    nrow = 4, ncol = 5, byrow = TRUE
)
bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

Q <- c(length(pir), length(pic))

# Create the fitBipartite objects
fitCol <- fitBipartiteSBMPop$new(A = bipartite_collection, Q = Q)


# Clustering by row and cols separately
fitCol.spectral_biclustering <- lapply(
    seq_along(fitCol$A),
    function(m) {
        spectral_biclustering(X = fitCol$A[[m]], K = fitCol$Q)
    }
)

fitCol.spectral_biclustering

# Co-Clustering

fitCol.spectral_coclustering <- lapply(
    seq_along(fitCol$A),
    function(m) {
        spectral_coclustering(X = fitCol$A[[m]], K = fitCol$Q)
    }
)

fitCol.spectral_coclustering


# # Add NAs to the collections
# colWithNA <- col
# colWithNA$A <- lapply(
#     seq_along(colWithNA$A),
#     function(m) {
#         colWithNA$A[[m]][sample(length(colWithNA$A[[m]]), numberOfNAsPerNetwork)] <- rep(NA, numberOfNAsPerNetwork)
#         colWithNA$A[[m]]
#     }
# )
# fitColWithNA <- fitBipartiteSBMPop$new(A = colWithNA$A, Q = Q)

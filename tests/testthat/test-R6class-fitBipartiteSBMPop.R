# Sourcing all necessary files
require("sbm")
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

set.seed(1234)
coupleOfNodesPerClasses <- c(10, 25)
Q <- c(3, 6)
numberOfNetworks <- 5
numberOfNAsPerNetwork <- 10

# Generate a collection
col <- generate_bipartite_collection(coupleOfNodesPerClass = coupleOfNodesPerClasses, coupleOfNumberOfClasses = Q, numberOfNetworks = numberOfNetworks)

# Add NAs to the collections
colWithNA <- col
colWithNA$A <- lapply(
    seq_along(colWithNA$A),
    function(m) {
        colWithNA$A[[m]][sample(length(colWithNA$A[[m]]), numberOfNAsPerNetwork)] <- rep(NA, numberOfNAsPerNetwork)
        colWithNA$A[[m]]
    }
)

# Create the fitBipartite objects
fitCol <- fitBipartiteSBMPop$new(A = col$A, Q = Q)
fitColWithNA <- fitBipartiteSBMPop$new(A = colWithNA$A, Q = Q)

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

eps <- 0.05
# Unipartite test of generation
# pi <- c(0.1, 0.6, 0.3)

# alphaCompartiment <- matrix(
#     c(
#         0.8, eps, 0.5,
#         eps, 0.5, eps,
#         eps, eps, .8
#     ),
#     nrow = 3
# )
# adj1 <- generate_unipartite_network(n = 50, pi = pi, alpha = alphaCompartiment)
# plotMyMatrix(adj1)
# sbm1 <- estimateSimpleSBM(adj1, estimOptions = list(plot=F, verbosity=1))
# plot(sbm1)

# Bipartite test of generation

nr <- 500
nc <- 2000

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

incid1 <- generate_bipartite_network(nr, nc, pir, pic, alpha)
lbm1 <- estimateBipartiteSBM(incid1)
plot(lbm1)

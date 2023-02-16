# Sourcing all necessary files
require("sbm")
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

set.seed(1234)

# TODO : essayer de débugger
# D'abord :
# Bien marquer les pir et pic pour ne pas compliquer la tâche
# Tirer au sort les alpha dans une uniforme
# alpha de dimensions 3x4
#
# Vérifier le clustering fournit par spectral_biclustering
# en utilisant aricode::ARI()
#
# Changer les prob utilisées pour permuter les lignes et colonnes


# Generate a bipartite collection fixed nr and nc

eps <- 0.05
M <- 5
nr <- 100
nc <- 250

pir <- rep(1/4, 4)
pic <- rep(1/5, 5)

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
fitCol <- fitBipartiteSBMPop$new(
    A = bipartite_collection,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    fit_opts = list(verbosity = 4)
)


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

sbm1 <- estimateBipartiteSBM(bipartite_collection[[1]])

fitCol$optimize()


# Add NAs to the collections
# bipartite_collection_NA <- bipartite_collection
# bipartite_collection_NA <- lapply(
#     seq_along(bipartite_collection_NA),
#     function(m) {
#         bipartite_collection_NA[[m]][sample(length(bipartite_collection_NA[[m]]), numberOfNAsPerNetwork)] <- rep(NA, numberOfNAsPerNetwork)
#         bipartite_collection_NA[[m]]
#     }
# )
# fitColWithNA <- fitBipartiteSBMPop$new(A = , Q = Q)

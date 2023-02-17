# Sourcing all necessary files
require("sbm", quietly = T)
require("ggpubr", quietly = T)
require("visNetwork", quietly = T)
require("igraph", quietly = T)
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")

set.seed(1234)

# TODO : essayer de débugger
# Changer les prob utilisées pour permuter les lignes et colonnes

# Generate a bipartite collection fixed nr and nc

eps <- 0.05
M <- 5
nr <- 100
nc <- 250

pir <- c(0.2, 0.8)
pic <- c(0.2, 0.3, 0.5)

Q <- c(length(pir), length(pic))

alpha <- matrix(runif(Q[1] * Q[2]), nrow = Q[1], ncol = Q[2])

ouput_file <- "output.txt"

cat(file = ouput_file, append = TRUE, "New run : \nParameters are : \n - M : ", M, "\n - Pi row : ", pir, "\n - Pi col : ", pic, "\n - alpha : ", alpha, "\n")


# alpha <- matrix(
#     c(
#         0.9, eps, eps, eps,
#         eps, 0.9, eps, eps,
#         eps, eps, 0.9, eps
#     ),
#     nrow = 3, ncol = 4, byrow = TRUE
# )

bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m){
    bipartite_collection[[m]]$incidence_matrix
})

## Init given with exact membership

Z <- lapply(seq.int(M), function(m){
    list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})

fitColExactMembership <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    Z = Z,
    free_density = FALSE,
    init_method = "given",
    fit_opts = list(verbosity = 4)
)

fitColExactMembership$optimize()

cat(file = ouput_file, append = TRUE, "Given exact clusters :\n - Real alpha : ", alpha, "\n - Given alpha : ", fitColExactMembership$alpha, "\n")

ARImethod <- function(fit_col){
    out <- list()
    out <- lapply(seq.int(fit_col$M),function(m){
        list(
            row = aricode::ARI(apply(fitColSpectral$tau[[m]][[1]], 1, which.max), Z[[m]][[1]]),
            col = aricode::ARI(apply(fitColSpectral$tau[[m]][[2]], 1, which.max), Z[[m]][[2]])
        )
    })
}

## Init with given LBM
sepLBM <- lapply(seq_along(bipartite_collection),
function(m) {
    estimateBipartiteSBM(bipartite_collection_incidence[[m]],
        estimOptions = list(
            verbosity = 1,
            plot = FALSE
        )
    )
})

# Extracting the list of the memberships
Z_LBM <- lapply(seq.int(M), function(m) {
    list(sepLBM[[m]]$memberships$row, sepLBM[[m]]$memberships$col)
})

LBM_reorganized_matrices <- lapply(seq.int(M), function(m){

})

# Create the fitBipartite object
fitColSepLBM <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    Z = Z_LBM,
    init_method = "given",
    fit_opts = list(verbosity = 4)
)


# Checking the clustering obtained
fitColSepLBM.resultARI <- as.data.frame(do.call("rbind", lapply(
    invisible(seq.int(M)),
    function(m) {
        ARIrow <- aricode::ARI(
            bipartite_collection[[m]]$row_clustering,
            sepLBM[[m]]$memberships$row
        )
        ARIcol <- aricode::ARI(
            bipartite_collection[[m]]$col_clustering,
            sepLBM[[m]]$memberships$col
        )
        list(step = m, ARIrow = ARIrow, ARIcol = ARIcol)
    }
)))

cat(file=ouput_file,append=TRUE,"SepLBM Clustering\n")
knitr::kable(fitColSepLBM.resultARI)

# SepLBM optimize
fitColSepLBM$optimize()
cat("SepLBM")
print(ARImethod(fitColSepLBM))

cat(file = ouput_file, append = TRUE, "SepLBM clustering :\n - Real alpha : ", alpha, "\n - SepLBM alpha : ", fitColSepLBM$alpha, "\n")

# Init with Spectral

# Create the fitBipartite object
fitColSpectral <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    fit_opts = list(verbosity = 4)
)

# Spectral Clustering by row and cols separately
fitColSpectral.spectralbiclust <- lapply(
    seq_along(fitColSpectral$A),
    function(m) {
        spectral_biclustering(X = fitColSpectral$A[[m]], K = fitColSpectral$Q)
    }
)

# Checking the clustering obtained
fitColSpectral.resultARI <- as.data.frame(do.call("rbind", lapply(
    invisible(seq.int(M)),
    function(m) {
        ARIrow <- aricode::ARI(
            bipartite_collection[[m]]$row_clustering,
            fitColSpectral.spectralbiclust[[m]]$row_clustering
        )
        ARIcol <- aricode::ARI(
            bipartite_collection[[m]]$col_clustering,
            fitColSpectral.spectralbiclust[[m]]$col_clustering
        )
        list(step = m, ARIrow = ARIrow, ARIcol = ARIcol)
    }
)))

cat(file=ouput_file,append=TRUE,"Spectral Bi-Clustering\n")
knitr::kable(fitColSpectral.resultARI)

# Spectral optimize
fitColSpectral$optimize()
cat("Spectral")
print(ARImethod(fitColSpectral))


cat(file=ouput_file,append=TRUE,"Spectral clustering :\n - Real alpha : ", alpha, "\n - Spectral alpha : ", fitColSpectral$alpha, "\n")

## Init with HCA
# Create the fitBipartite object
fitColHCA <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q, free_mixture = FALSE,
    free_density = FALSE,
    init_method = "hca",
    fit_opts = list(verbosity = 4)
)

# Hierarchical Clustering by row and cols separately
fitColHCA.hierarchicalbiclust <- lapply(
    seq_along(fitColHCA$A),
    function(m) {
        bipartite_hierarchic_clustering(X = fitColHCA$A[[m]], K = fitColHCA$Q)
    }
)

# Checking the clustering obtained
hierarchicalbiclust.resultARI <- as.data.frame(do.call("rbind", lapply(
    invisible(seq.int(M)),
    function(m) {
        ARIrow <- aricode::ARI(
            bipartite_collection[[m]]$row_clustering,
            fitColHCA.hierarchicalbiclust[[m]]$row_clustering
        )
        ARIcol <- aricode::ARI(
            bipartite_collection[[m]]$col_clustering,
            fitColHCA.hierarchicalbiclust[[m]]$col_clustering
        )
        list(step = m, ARIrow = ARIrow, ARIcol = ARIcol)
    }
)))

cat(file=ouput_file,append=TRUE,"Hierarchical Bi-Clustering\n")
knitr::kable(hierarchicalbiclust.resultARI)

fitColHCA$optimize()
cat("HCA")
print(ARImethod(fitColSpectral))


cat(file=ouput_file,append=TRUE,"HCA clustering :\n - Real alpha : ", alpha, "\n - HCA alpha : ", fitColHCA$alpha, "\n\n")

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

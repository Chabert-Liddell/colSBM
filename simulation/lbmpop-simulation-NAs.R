# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

devtools:::load_all(path = "R/")

set.seed(1234)

verbose <- TRUE
test_alea <- TRUE

eps <- 0.05
M <- 3
nr <- 100
nc <- 250

pir <- c(0.2, 0.8)
pic <- c(0.2, 0.3, 0.5)

Q <- c(length(pir), length(pic))

alpha <- matrix(
  c(
    0.9, eps, 0.5,
    eps, 0.8, 0.2
  ),
  nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M, return_memberships = TRUE)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
  bipartite_collection[[m]]$incidence_matrix
})
NAs_index <- sample(seq_len(length(bipartite_collection_incidence[[1]])), floor(0.5 * length(bipartite_collection_incidence[[1]])))
real_val_NAs <- bipartite_collection_incidence[[1]][NAs_index]
bipartite_collection_incidence[[1]][NAs_index] <- NA
NAs_coordinates <- which(is.na(bipartite_collection_incidence[[1]]), arr.ind = TRUE)
x_NAs <- sort(unique(NAs_coordinates[, 1]))
y_NAs <- sort(unique(NAs_coordinates[, 2]))
## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
  list(bipartite_collection[[m]]$row_blockmemberships, bipartite_collection[[m]]$col_blockmemberships)
})
tic()
mybisbmpop <- estimate_colBiSBM(
  netlist = bipartite_collection_incidence, colsbm_model = "iid",
  global_opts = list(
    # parallelization_vector = c(F,F),
    nb_cores = 6, verbosity = 4
  )
)
toc()
# choosed_bisbmpop <- estimate_colBiSBM(
#     netlist = bipartite_collection_incidence,
#     colsbm_model = "iid",
#     global_opts = list(nb_cores = 3)
# )

ari_sums <- sapply(
  seq_along(mybisbmpop$best_fit$Z),
  function(m) {
    sum(c(
      aricode::ARI(
        Z[[m]][[1]],
        mybisbmpop$best_fit$Z[[m]][[1]]
      ),
      aricode::ARI(
        Z[[m]][[2]],
        mybisbmpop$best_fit$Z[[m]][[2]]
      )
    ))
  }
)

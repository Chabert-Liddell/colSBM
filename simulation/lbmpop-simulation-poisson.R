# Sourcing all necessary files
require("sbm", quietly = T)
require("dplyr", quietly = T)
require("tictoc", quietly = T)
require("ggplot2", quietly = T)

devtools::load_all(path = "R/")

set.seed(1234)

verbose <- TRUE
test_alea <- TRUE

eps <- 0.05
M <- 3
nr <- 100
nc <- 250

pic <- c(0.1, 0.2, 0.3, 0.4)
pir <- c(0.1, 0.2, 0.7)

Q <- c(length(pir), length(pic))

alpha <- matrix(
  c(
    1, 0, 0, 5,
    4, 100, 25, 45,
    0, 55, 0, 60
  ),
  nrow = Q[1], ncol = Q[2], byrow = TRUE
)

bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M, distribution = "poisson", return_memberships = TRUE)

# This is a list of the M incidence matrices
bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
  bipartite_collection[[m]]$incidence_matrix
})

bipartite_collection_incidence_binary <- lapply(seq_along(bipartite_collection_incidence), function(m) {
  1 * (bipartite_collection_incidence[[m]] > 0)
})

## Init given with exact membership

Z <- lapply(seq.int(M), function(m) {
  list(bipartite_collection[[m]]$row_clustering, bipartite_collection[[m]]$col_clustering)
})
tic()
mybisbmpop <- estimate_colBiSBM(
  netlist = bipartite_collection_incidence_binary, colsbm_model = "iid",
  nb_run = 3,
  distribution = "poisson",
  global_opts = list(
    parallelization_vector = c(T, T),
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

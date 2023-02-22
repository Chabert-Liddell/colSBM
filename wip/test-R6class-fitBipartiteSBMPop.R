require("testthat")
require("tictoc")
source("R/R6class-fitBipartiteSBMPop.R")
source("R/utils.R")

# Parallelize
library(doParallel)

# Setup backend to use many processors
totalCores <- detectCores()

# Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1] - 1)
registerDoParallel(cluster)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

set.seed(1234)

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


iter_max <- 100

vanilla.ICL <- rep(-Inf, iter_max)

minibatch.ICL <- rep(-Inf, iter_max)
tic("All loops")
for (iter in 1:iter_max) {
  set
  tic("Vanilla")
  fitCol <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q,
    free_mixture = FALSE,
    free_density = FALSE,
    init_method = "spectral",
    fit_opts = list(verbosity = 0, minibatch = FALSE)
  )
  fitCol$optimize()
  toc(log = TRUE, quiet = TRUE)
  vanilla.ICL[iter] <- fitCol$MAP$ICL
}

vanilla.lst <- tic.log(format = FALSE)

vanilla.timings <- unlist(lapply(vanilla.lst, function(x) x$toc - x$tic))
mean(vanilla.timings)

for (iter in 1:iter_max) {
  tic("Minibatch")
  fitColMinibatch <- fitBipartiteSBMPop$new(
    A = bipartite_collection_incidence,
    Q = Q,
    free_mixture = FALSE,
    free_density = FALSE,
    init_method = "spectral",
    fit_opts = list(verbosity = 0, minibatch = TRUE)
  )
  fitColMinibatch$optimize()
  toc(log = TRUE, quiet = TRUE)
  minibatch.ICL[iter] <- fitColMinibatch$MAP$ICL
}

minibatch.lst <- tic.log(format = FALSE)
minibatch.timings <- unlist(lapply(minibatch.lst, function(x) x$toc - x$tic))
mean(minibatch.timings)
toc()
stopCluster(cluster)

vanilla.ICL
minibatch.ICL

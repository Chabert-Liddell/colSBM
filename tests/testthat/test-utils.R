if (!require("testthat")) install.packages("testthat")
library("testthat")
source("R/utils.R")
# Testing the generate_bipartite_collection
set.seed(1234)

eps <- 0.05

pir <- c(0.1, 0.5, 0.2, 0.2)
pic <- c(1 / 3, 1 / 3, 1 / 3)

alpha <- matrix(c(
    0.7, 0.6, 0.6,
    0.6, 0.5, eps,
    0.6, eps, eps,
    eps, eps, 0.9
), nrow = 4, byrow = TRUE)

M <- 5

# With one nr and one nc
nr <- 120
nc <- 240

firstCol <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("nc and nr are extended to match M", {
    expect_equal(length(firstCol), M)
    expect_equal(unique(lapply(firstCol, nrow)), list(nr))
    expect_equal(unique(lapply(firstCol, ncol)), list(nc))
})

# With increasing nr and nc
nr <- c(5, 10, 20, 40, 80)
nc <- 2 * nr

increasingCol <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("Generate bipartite collection with vectors for nc and nr", {
    expect_equal(length(increasingCol), M)
    expect_equal(lapply(increasingCol, nrow), as.list(nr))
    expect_equal(lapply(increasingCol, ncol), as.list(nc))
})

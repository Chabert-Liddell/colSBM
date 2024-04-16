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

firstCol <- colSBM::generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("nc and nr are extended to match M", {
    expect_equal(length(firstCol), M)
    expect_equal(unique(lapply(firstCol, nrow)), list(nr))
    expect_equal(unique(lapply(firstCol, ncol)), list(nc))
})

# With increasing nr and nc
nr <- c(5, 10, 20, 40, 80)
nc <- 2 * nr

increasingCol <- colSBM::generate_bipartite_collection(nr, nc, pir, pic, alpha, M)

test_that("Generate bipartite collection with vectors for nc and nr", {
    expect_equal(length(increasingCol), M)
    expect_equal(lapply(increasingCol, nrow), as.list(nr))
    expect_equal(lapply(increasingCol, ncol), as.list(nc))
    #  Wrong M
    expect_error(colSBM::generate_bipartite_collection(
        nr, nc, pir, pic, alpha,
        M + 1
    ))
    expect_error(colSBM::generate_bipartite_collection(
        c(nr, 1), nc, pir, pic, alpha,
        M
    ))
    expect_error(colSBM::generate_bipartite_collection(
        nr, c(nc, 1), pir, pic, alpha,
        M
    ))
})

test_that("Testing the different models", {
    expect_error(
        colSBM::generate_bipartite_collection(
            nr, nc, pir, pic, alpha,
            M,
            model = "iidi"
        )
    )
    expect_no_error(
        colSBM::generate_bipartite_collection(
            nr, nc, pir, pic, alpha,
            M,
            model = "iid"
        )
    )
    expect_no_error(
        colSBM::generate_bipartite_collection(
            nr, nc, pir, pic, alpha,
            M,
            model = "pi"
        )
    )
    expect_no_error(
        colSBM::generate_bipartite_collection(
            nr, nc, pir, pic, alpha,
            M,
            model = "rho"
        )
    )
    expect_no_error(
        colSBM::generate_bipartite_collection(
            nr, nc, pir, pic, alpha,
            M,
            model = "pirho"
        )
    )
})

test_that("Base case spectral biclustering", {
    X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
    expect_equal(
        colSBM:::spectral_biclustering(X = X, K = c(1, 1)),
        list(
            row_clustering = c(1, 1),
            col_clustering = c(1, 1)
        )
    )
})

test_that("Spectral clustering too many clusters works", {
    X <- generate_unipartite_network(n = 10, pi = 1, alpha = 0.9)
    expect_message(
        colSBM:::spectral_clustering(X = X, K = 12),
    )
})

test_that("Base case spectral biclustering", {
    X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
    expect_equal(
        colSBM:::spectral_biclustering(X = X, K = c(1, 1)),
        list(
            row_clustering = c(1, 1),
            col_clustering = c(1, 1)
        )
    )
})

test_that("Base case hierarchical biclustering", {
    X <- matrix(c(1, 0, 0, 1), byrow = TRUE, nrow = 2)
    expect_equal(
        colSBM:::bipartite_hierarchic_clustering(X = X, K = c(1, 1)),
        list(
            row_clustering = c(1, 1),
            col_clustering = c(1, 1)
        )
    )
})

test_that("Hierarchical biclustering works", {
    X <- generate_bipartite_network(
        nr = 10, 10, pi = c(0.4, 0.6), rho = 1,
        alpha = matrix(c(0.9, 0.1), byrow = TRUE, nrow = 2)
    )
    expect_no_error(
        colSBM:::bipartite_hierarchic_clustering(X = X, K = c(2, 1)),
    )
})

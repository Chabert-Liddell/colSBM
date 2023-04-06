# iid model test
set.seed(1234)

iid_expected_BICL <- -30533.22

iid_eps <- 0.05
iid_M <- 3
iid_nr <- 100
iid_nc <- 250
iid_pir <- c(0.2, 0.8)
iid_pic <- c(0.2, 0.8)
iid_Q <- c(length(iid_pir), length(iid_pic))
iid_alpha <- matrix(
    c(
        0.9, iid_eps,
        iid_eps, 0.8
    ),
    nrow = iid_Q[1], ncol = iid_Q[2], byrow = TRUE
)

iid_bipartite_collection <-
    colSBM::generate_bipartite_collection(iid_nr, iid_nc, iid_pir, iid_pic, iid_alpha, iid_M)

# This is a list of the M incidence matrices
iid_bipartite_collection_incidence <- lapply(seq.int(iid_M), function(m) {
    iid_bipartite_collection[[m]]$incidence_matrix
})

iid_Z <- lapply(seq.int(iid_M), function(m) {
    list(
        iid_bipartite_collection[[m]]$row_clustering,
        iid_bipartite_collection[[m]]$col_clustering
    )
})

iid_choosed_bisbmpop <- colSBM::estimate_colBiSBM(
    netlist = iid_bipartite_collection_incidence,
    colsbm_model = "iid",
    global_opts = list(nb_cores = 3, plot_details = 0)
)

iid_ari_sums <- sapply(
    seq_along(iid_choosed_bisbmpop$best_fit$Z),
    function(m) {
        sum(c(
            aricode::ARI(
                iid_Z[[m]][[1]],
                iid_choosed_bisbmpop$best_fit$Z[[m]][[1]]
            ),
            aricode::ARI(
                iid_Z[[m]][[2]],
                iid_choosed_bisbmpop$best_fit$Z[[m]][[2]]
            )
        ))
    }
)

test_that("iid model returns correct object", {
    expect_equal(class(iid_choosed_bisbmpop), c("bisbmpop", "R6"))
})

test_that("iid model returns correct values on simulation", {
    expect_equal(iid_ari_sums, rep(2,iid_M))
    expect_equal(iid_choosed_bisbmpop$best_fit$Q, iid_Q)
    expect_equal(round(iid_choosed_bisbmpop$best_fit$BICL,2), iid_expected_BICL)
    #TODO fix permutation of alpha : expect_equal(round(iid_choosed_bisbmpop$best_fit$alpha, 2), iid_alpha)
})

require("ggplot2")
require("tictoc")

devtools::load_all("R/")

# Generation of conditions
if (!exists("model_to_test")) {
    model_to_test <- "iid"
}

if (!exists("repetitions")) {
    repetitions <- seq.int(3)
}

nr <- 75
nc <- 75

pi <- matrix(c(0.2, 0.3, 0.5), nrow = 1, byrow = TRUE)
rho <- matrix(rep(1 / 3, 3), nrow = 1, byrow = TRUE)

if (!exists("arg")) {
    arg <- commandArgs(trailingOnly = TRUE)
}

if (identical(arg, character(0))) {
    epsilons <- c(0.1, 0.4)
} else {
    epsilons <- as.numeric(arg)
}

conditions <- tidyr::crossing(epsilons, pi, rho, repetitions)

results <- bettermc::mclapply(seq_len(nrow(conditions)), function(s) {
    eps <- conditions[s, ]$epsilons
    current_pi <- conditions[s, ]$pi
    current_rho <- conditions[s, ]$rho

    alpha_assortative <- matrix(0.3, nrow = 3, ncol = 3) +
        matrix(
            c(
                eps, -0.5 * eps, -0.5 * eps,
                -0.5 * eps, eps, -0.5 * eps,
                -0.5 * eps, -0.5 * eps, eps
            ),
            nrow = 3, byrow = TRUE
        )

    alpha_core_periphery <- matrix(0.3, nrow = 3, ncol = 3) +
        matrix(
            c(
                1.5 * eps, eps, 0.5 * eps,
                eps, 0.5 * eps, 0,
                0.5 * eps, 0, -0.5 * eps
            ),
            nrow = 3, byrow = TRUE
        )

    alpha_disassortative <- matrix(0.3, nrow = 3, ncol = 3) +
        matrix(
            c(
                -0.5 * eps, eps, eps,
                eps, -0.5 * eps, eps,
                eps, eps, -0.5 * eps
            ),
            nrow = 3, byrow = TRUE
        )

    assortative_collection <- generate_bipartite_collection(
        nr, nc,
        current_pi, current_rho,
        alpha_assortative, 3,
        model = model_to_test
    )

    assortative_incidence <- lapply(
        seq_along(assortative_collection),
        function(m) {
            return(assortative_collection[[m]]$incidence_matrix)
        }
    )

    assortative_row_clustering <- lapply(
        seq_along(assortative_collection),
        function(m) {
            return(assortative_collection[[m]]$row_clustering)
        }
    )

    assortative_col_clustering <- lapply(
        seq_along(assortative_collection),
        function(m) {
            return(assortative_collection[[m]]$row_clustering)
        }
    )

    core_periphery_collection <- generate_bipartite_collection(
        nr, nc,
        current_pi, current_rho,
        alpha_core_periphery, 3,
        model = model_to_test
    )

    core_periphery_incidence <- lapply(
        seq_along(core_periphery_collection),
        function(m) {
            return(core_periphery_collection[[m]]$incidence_matrix)
        }
    )

    core_periphery_row_clustering <- lapply(
        seq_along(core_periphery_collection),
        function(m) {
            return(core_periphery_collection[[m]]$row_clustering)
        }
    )

    core_periphery_col_clustering <- lapply(
        seq_along(core_periphery_collection),
        function(m) {
            return(core_periphery_collection[[m]]$row_clustering)
        }
    )

    disassortative_collection <- generate_bipartite_collection(
        nr, nc,
        current_pi, current_rho,
        alpha_disassortative, 3,
        model = model_to_test
    )

    disassortative_incidence <- lapply(
        seq_along(disassortative_collection),
        function(m) {
            return(disassortative_collection[[m]]$incidence_matrix)
        }
    )

    disassortative_row_clustering <- lapply(
        seq_along(disassortative_collection),
        function(m) {
            return(disassortative_collection[[m]]$row_clustering)
        }
    )

    disassortative_col_clustering <- lapply(
        seq_along(disassortative_collection),
        function(m) {
            return(disassortative_collection[[m]]$row_clustering)
        }
    )

    real_row_clustering <- append(
        append(
            assortative_row_clustering,
            core_periphery_row_clustering
        ),
        disassortative_row_clustering
    )

    real_col_clustering <- append(
        append(
            assortative_col_clustering,
            core_periphery_col_clustering
        ),
        disassortative_col_clustering
    )

    incidence_matrices <- append(
        append(
            assortative_incidence,
            core_periphery_incidence
        ),
        disassortative_incidence
    )

    netids <- rep(c("as", "cp", "dis"), each = 3)

    tic()
    list_collection <- clusterize_bipartite_networks(
        netlist = incidence_matrices,
        net_id = netids,
        colsbm_model = model_to_test,
        global_opts = list(
            nb_cores = parallel::detectCores() - 1, verbosity = 1,
            plot_details = 0
        ),
        silent_parallelization = TRUE
    )
    toc()
    return(
        list(
        epsilon = eps,
        repetition = repetition,
        list_of_clusterings = list_collection,
        real_block_memberships = list(
            row = real_row_clustering,
            col = real_col_clustering
        )
    ))
},
mc.cores = parallel::detectCores() - 1,
mc.progress = TRUE,
mc.retry = -1
)

saveRDS(results, file = paste0(
    "simulation/data/",
    "simulated_collection_clustering_",
    model_to_test, "_",
    format(Sys.time(), "%d-%m-%y-%X"),
    ".Rds"
))

require("bettermc")
require("gtools")
require("tictoc")
devtools::load_all("R/")

# Network param
nr <- 120
nc <- 120
M <- 2

# Changing parameters
base_alpha <- matrix(0.25, nrow = 4, ncol = 4)
epsilon_alpha <- seq(from = 0.0, to = 0.24, by = 0.03)

base_pi1 <- c(0.2, 0.4, 0.4, 0)
rho1 <- rep(0.25, 4)

pi2 <- rep(0.25, 4)
base_rho2 <- c(0, 1/3, 1/3, 1/3)

pi1 <- matrix(unlist(combinat::permn(base_pi1)), byrow = TRUE, ncol = 4)
pi1 <- pi1[!duplicated(pi1), ][1:4, ]

rho2 <- matrix(unlist(combinat::permn(base_rho2)), byrow = TRUE, ncol = 4)
rho2 <- rho2[!duplicated(rho2),]

repetition <- seq.int(3)

conditions <- tidyr::crossing(epsilon_alpha, pi1, rho2, repetition)

# Filter conditions to prevent the same blocks from being empty
conditions <- conditions[
    !apply(conditions$pi1[, 1:4] == 0 & conditions$rho2[, 1:4] == 0, 
    1, any),
]

# To speed up computations and debug adding an argument based selection
arg <- commandArgs(trailingOnly = TRUE)
if (identical(arg, character(0))) {
    cat(
        "\nNo arguments provided,",
        "assuming you want to go over all the conditions."
    )
    arg <- c(1, nrow(conditions))
} else {
    arg <- as.numeric(arg)
}
if (arg[1] < 1 | arg[1] > nrow(conditions)) {
    warning(paste("Arg 1 was invalid, set to 1."))
    arg[1] <- 1
}
if (arg[2] > nrow(conditions) | arg[2] < 1) {
    warning(paste("Arg 2 was invalid, set to", nrow(conditions)))
    arg[2] <- nrow(conditions)
}

choosed_conditions <- seq.int(from = arg[1], to = arg[2])

conditions <- conditions[choosed_conditions,]
tic()
results <- bettermc::mclapply(seq_len(nrow(conditions)), function(s) {
    ea <- conditions[s,]$epsilon_alpha
    current_pi1 <- conditions[s, ]$pi1
    current_rho2 <- conditions[s,]$rho2
    
    current_alpha <- base_alpha + matrix(c(
                            3 * ea, 2 * ea, ea, -ea,
                            2 * ea, 2 * ea, - ea, ea,
                            ea, - ea, ea, 2 * ea,
                            - ea, ea, 2 * ea, 0),
        byrow = TRUE, nrow = 4, ncol = 4
    )

    # Compute supports
    Cpi1 <- matrix(c(current_pi1, pi2), byrow = TRUE, nrow = M) > 0
    Cpi2 <- matrix(c(rho1, current_rho2), byrow = TRUE, nrow = M) > 0

    netlist_generated <- list(
        generate_bipartite_network(
            nr, nc, conditions[s, ]$pi1, rho1,
            current_alpha
        ),
        generate_bipartite_network(
            nr, nc, pi2, conditions[s, ]$rho2,
            current_alpha
        )
    )
    netlist <- lapply(seq_along(netlist_generated), function(m) {
        return(netlist_generated[[m]]$incidence_matrix)
    })

    row_clusterings <- lapply(seq_along(netlist_generated), function(m) {
        return(netlist_generated[[m]]$row_clustering)
    })

    col_clusterings <- lapply(seq_along(netlist_generated), function(m) {
        return(netlist_generated[[m]]$col_clustering)
    })

    full_row_clustering <- as.vector(sapply(
        seq.int(M),
        function(m) row_clusterings[[m]]
    ))

    full_col_clustering <- as.vector(sapply(
        seq.int(M),
        function(m) col_clusterings[[m]]
    ))

    fitted_bisbmpop_iid <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "iid",
	nb_run = 1,
        silent_parallelization = TRUE,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_pi <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "pi",
	nb_run = 1,
        silent_parallelization = TRUE,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_rho <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "rho",
	nb_run = 1,
        silent_parallelization = TRUE,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_pirho <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "pirho",
	nb_run = 1,
        silent_parallelization = TRUE,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    # BICLs
    sep_BICL <- sum(fitted_bisbmpop_iid$sep_BiSBM$BICL)
    iid_BICL <- fitted_bisbmpop_iid$best_fit$BICL
    pi_BICL <- fitted_bisbmpop_pi$best_fit$BICL
    rho_BICL <- fitted_bisbmpop_rho$best_fit$BICL
    pirho_BICL <- fitted_bisbmpop_pirho$best_fit$BICL
    BICLs <- c(sep_BICL, iid_BICL, pi_BICL, rho_BICL, pirho_BICL)

    # ARIs
    compute_mean_ARI <- function(model) {
        # We compute the mean amongst the two networks and return values for
        # rows and columns in a vector
        # sapply ives a matrix with in row the axis ARIs
        # and in columns the networks
        #    1     2
        # ax row1  row2
        # ay col1  col2
        rowMeans(sapply(seq.int(model$M), function(m) {
            c(
                aricode::ARI(model$Z[[m]][[1]], row_clusterings[[m]]),
                aricode::ARI(model$Z[[m]][[2]], col_clusterings[[m]])
            )
        }))
    }
    sep_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_iid$sep_BiSBM)
    iid_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_iid$best_fit)
    pi_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_pi$best_fit)
    rho_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_rho$best_fit)
    pirho_mean_ARIs <- compute_mean_ARI(fitted_bisbmpop_pirho$best_fit)

    compute_double_ARI <- function(model) {
        model_row_Z <- as.vector(sapply(
            seq.int(model$M),
            function(m) model$Z[[m]][[1]]
        ))

        model_col_Z <- as.vector(sapply(
            seq.int(model$M),
            function(m) model$Z[[m]][[2]]
        ))

        return(list(
            aricode::ARI(model_row_Z, full_row_clustering),
            aricode::ARI(model_col_Z, full_col_clustering)
        ))
    }

    sep_double_ARIs <- compute_double_ARI(fitted_bisbmpop_iid$sep_BiSBM)
    iid_double_ARIs <- compute_double_ARI(fitted_bisbmpop_iid$best_fit)
    pi_double_ARIs <- compute_double_ARI(fitted_bisbmpop_pi$best_fit)
    rho_double_ARIs <- compute_double_ARI(fitted_bisbmpop_rho$best_fit)
    pirho_double_ARIs <- compute_double_ARI(fitted_bisbmpop_pirho$best_fit)

    data_frame_output <- data.frame(
        # The conditions
        epsilon_alpha = ea,
        pi1 = current_pi1,
        rho2 = current_rho2,
        repetition = as.numeric(conditions[s, 4]),
        # The results
        ## sep
        sep_BICL = sep_BICL,
        sep_mean_row_ARI = sep_mean_ARIs[1],
        sep_mean_col_ARI = sep_mean_ARIs[2],
        sep_double_row_ARI = sep_double_ARIs[[1]],
        sep_double_col_ARI = sep_double_ARIs[[2]],

        ## iid
        iid_BICL = iid_BICL,
        iid_mean_row_ARI = iid_mean_ARIs[1],
        iid_mean_col_ARI = iid_mean_ARIs[2],
        iid_double_row_ARI = iid_double_ARIs[[1]],
        iid_double_col_ARI = iid_double_ARIs[[2]],
        iid_Q1 = fitted_bisbmpop_iid$best_fit$Q[1],
        iid_Q2 = fitted_bisbmpop_iid$best_fit$Q[2],

        ## pi
        pi_BICL = pi_BICL,
        pi_mean_row_ARI = pi_mean_ARIs[1],
        pi_mean_col_ARI = pi_mean_ARIs[2],
        pi_double_row_ARI = pi_double_ARIs[[1]],
        pi_double_col_ARI = pi_double_ARIs[[2]],
        pi_Q1 = fitted_bisbmpop_pi$best_fit$Q[1],
        pi_Q2 = fitted_bisbmpop_pi$best_fit$Q[2],

        ## pi
        rho_BICL = rho_BICL,
        rho_mean_row_ARI = rho_mean_ARIs[1],
        rho_mean_col_ARI = rho_mean_ARIs[2],
        rho_double_row_ARI = rho_double_ARIs[[1]],
        rho_double_col_ARI = rho_double_ARIs[[2]],
        rho_Q1 = fitted_bisbmpop_rho$best_fit$Q[1],
        rho_Q2 = fitted_bisbmpop_rho$best_fit$Q[2],

        ## pirho
        pirho_BICL = pirho_BICL,
        pirho_mean_row_ARI = pirho_mean_ARIs[1],
        pirho_mean_col_ARI = pirho_mean_ARIs[2],
        pirho_double_row_ARI = pirho_double_ARIs[[1]],
        pirho_double_col_ARI = pirho_double_ARIs[[2]],
        pirho_Q1 = fitted_bisbmpop_pirho$best_fit$Q[1],
        pirho_Q2 = fitted_bisbmpop_pirho$best_fit$Q[2]
    )

    return(data_frame_output)
},
mc.cores = parallel::detectCores() - 1,
mc.progress = TRUE,
mc.retry = -1
)
toc()
full_data_frame <- do.call(rbind, results)

saveRDS(full_data_frame,
    file = paste0(
        "./simulation/data/inference_testing_",
        Sys.time(), "_", toString(choosed_conditions), ".Rds"
    )
)

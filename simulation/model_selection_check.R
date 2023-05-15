require("bettermc")
require("gtools")
require("tictoc")
devtools::load_all("R/")

# Network param
nr <- 90
nc <- 90

# Changing parameters
epsilons_pi <- seq(from = 0.0, to = 0.28, by = 0.035)
epsilons_rho <- seq(from = 0.0, to = 0.28, by = 0.035)
pi1 <- matrix(rep(1 / 3, 3), nrow = 1)
rho1 <- matrix(rep(1 / 3, 3), nrow = 1)
ea <- 0.16
alpha <- 0.25 + matrix(
    c(
        3 * ea, 2 * ea, ea,
        2 * ea, 2 * ea, -ea,
        ea,     -ea,    ea
    ),
    byrow = TRUE, nrow = 3, ncol = 3
)


prob_order <- seq(1,3)
prob_order <- t(sapply(combinat::permn(prob_order), function(v) v))


repetitions <- seq.int(3)

conditions <- tidyr::crossing(
    epsilon_pi = epsilons_pi,
    epsilon_rho = epsilons_rho,
    pi2_order = prob_order,
    rho2_order = prob_order,
    repetition = repetitions
)

# To speed up computations and debug adding an argument based selection
if (!exists("arg")) {
    arg <- commandArgs(trailingOnly = TRUE)
}
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

conditions <- conditions[choosed_conditions, ]
tic()
results <- bettermc::mclapply(seq_len(nrow(conditions)), function(s) {
    epsilon_pi <- conditions[s,]$epsilon_pi
    epsilon_rho <- conditions[s, ]$epsilon_rho

    # Computing the vector with the epsilons
    current_pi2 <- c(
        1 / 3 - epsilon_pi,
        1 / 3,
        1 / 3 + epsilon_pi
    )
    current_rho2 <- c(
        1 / 3 - epsilon_rho,
        1 / 3,
        1 / 3 + epsilon_rho
    )

    # Permutating the vectors
    current_pi2 <- current_pi2[conditions[s, ]$pi2_order]
    current_rho2 <- current_rho2[conditions[s, ]$rho2_order]

    netlist_generated <- list(
        generate_bipartite_network(
            nr, nc, pi1, rho1,
            alpha
        ),
        generate_bipartite_network(
            nr, nc, current_pi2, current_rho2,
            alpha
        )
    )

    # Extracting the incidence matrices
    netlist <- lapply(seq_along(netlist_generated), function(m) {
        return(netlist_generated[[m]]$incidence_matrix)
    })

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

    data_frame_output <- data.frame(
        # The conditions
        epsilon_pi = epsilon_pi,
        epsilon_rho = epsilon_rho,
        pi2 = matrix(current_pi2, nrow = 1),
        rho2 = matrix(current_rho2, nrow = 1),
        repetition = as.numeric(conditions[s,]$repetition),

        # The results
        ## sep
        sep_BICL = sep_BICL,

        ## iid
        iid_BICL = iid_BICL,
        iid_Q1 = fitted_bisbmpop_iid$best_fit$Q[1],
        iid_Q2 = fitted_bisbmpop_iid$best_fit$Q[2],

        ## pi
        pi_BICL = pi_BICL,
        pi_Q1 = fitted_bisbmpop_pi$best_fit$Q[1],
        pi_Q2 = fitted_bisbmpop_pi$best_fit$Q[2],

        ## pi
        rho_BICL = rho_BICL,
        rho_Q1 = fitted_bisbmpop_rho$best_fit$Q[1],
        rho_Q2 = fitted_bisbmpop_rho$best_fit$Q[2],

        ## pirho
        pirho_BICL = pirho_BICL,
        pirho_Q1 = fitted_bisbmpop_pirho$best_fit$Q[1],
        pirho_Q2 = fitted_bisbmpop_pirho$best_fit$Q[2],

        # Preferred model
        preferred_model = c("sep", "iid", "pi", "rho", "pirho")[which.max(BICLs)]
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
        "./simulation/data/model_selection_check_", max(repetitions), "_rep_",
        toString(sprintf("%04d", arg)), "_", 
        format(Sys.time(), "%d-%m-%y_%H-%M"),
        ".Rds"
    )
)

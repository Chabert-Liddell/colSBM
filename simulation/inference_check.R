require("bettermc")
require("gtools")
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
pi1 <- pi1[!duplicated(pi1),]

rho2 <- matrix(unlist(combinat::permn(base_rho2)), byrow = TRUE, ncol = 4)
rho2 <- rho2[!duplicated(rho2),]

model <- c("sep","iid", "pi", "rho", "pirho")

repetition <- seq.int(3)

conditions <- tidyr::crossing(epsilon_alpha, pi1, rho2, repetition)

# Filter conditions to prevent the same blocks from being empty
conditions <- conditions[
    !apply(conditions$pi1[, 1:4] == 0 & conditions$rho2[, 1:4] == 0, 
    1, any),
]

results <- bettermc::mclapply(seq_len(nrow(conditions)), function(c) {
    ea <- conditions[c,]$epsilon_alpha
    pi1 <- conditions[c, ]$pi1
    rho2 <- conditions[c,]$rho2
    
    current_alpha <- base_alpha + matrix(c(
                            3 * ea, 2 * ea, ea, -ea,
                            2 * ea, 2 * ea, - ea, ea,
                            ea, - ea, ea, 2 * ea,
                            - ea, ea, 2 * ea, 0),
        byrow = TRUE, nrow = 4, ncol = 4
    )

    # Compute supports
    Cpi1 <- matrix(c(pi1, pi2), byrow = TRUE, nrow = M) > 0
    Cpi2 <- matrix(c(rho1, rho2), byrow = TRUE, nrow = M) > 0

    netlist_generated <- list(
        generate_bipartite_network(
            nr, nc, conditions[c, ]$pi1, rho1,
            current_alpha
        ),
        generate_bipartite_network(
            nr, nc, pi2, conditions[c, ]$rho2,
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

    fitted_bisbmpop_iid <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "iid",
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_pi <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "pi",
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_rho <- estimate_colBiSBM(
        netlist = "rho",
        colsbm_model = conditions[c, ]$model,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = parallel::detectCores() - 1
        )
    )

    fitted_bisbmpop_pirho <- estimate_colBiSBM(
        netlist = netlist,
        colsbm_model = "pirho",
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
    compute_ARI <- function(model, Z) {
        # We compute the mean amongst the two networks and return values for
        # rows and columns in a vector
        rowMeans(sapply(seq.int(model$M), function(m) {
            c(
                aricode::ARI(model$Z[[m]][[1]], row_clusterings[[m]]),
                aricode::ARI(model$Z[[m]][[2]], col_clusterings[[m]])
            )
        }))
    }
    sep_ARIs <- compute_ARI(fitted_bisbmpop_iid$sep_BiSBM)
    iid_ARIs <- compute_ARI(fitted_bisbmpop_iid$best_fit)
    pi_ARIs <- compute_ARI(fitted_bisbmpop_pi$best_fit)
    rho_ARIs <- compute_ARI(fitted_bisbmpop_rho$best_fit)
    pirho_ARIs <- compute_ARI(fitted_bisbmpop_pirho$best_fit)
    row_ARIs <- c(
        sep_ARIs[1], iid_ARIs[1], pi_ARIs[1], 
        rho_ARIs[1], pirho_ARIs[1]
    )

    col_ARIs <- c(
        sep_ARIs[2], iid_ARIs[2], pi_ARIs[2],
        rho_ARIs[2], pirho_ARIs[2]
    )

    # Recovered blocks
    ## Less
    pi_Q1_less <- fitted_bisbmpop_pi$Q[1] < 4
    pi_Q2_less <- fitted_bisbmpop_pi$Q[2] < 4

    rho_Q1_less <- fitted_bisbmpop_rho$Q[1] < 4
    rho_Q2_less <- fitted_bisbmpop_rho$Q[2] < 4

    pirho_Q1_less <- fitted_bisbmpop_pirho$Q[1] < 4
    pirho_Q2_less <- fitted_bisbmpop_pirho$Q[2] < 4

    Q1_less <- c(NA, NA, pi_Q1_less, rho_Q1_less, pirho_Q1_less)
    Q2_less <- c(NA, NA, pi_Q2_less, rho_Q2_less, pirho_Q2_less)

    ## Greater
    pi_Q1_great <- fitted_bisbmpop_pi$Q[1] > 4
    pi_Q2_great <- fitted_bisbmpop_pi$Q[2] > 4

    rho_Q1_great <- fitted_bisbmpop_rho$Q[1] > 4
    rho_Q2_great <- fitted_bisbmpop_rho$Q[2] > 4

    pirho_Q1_great <- fitted_bisbmpop_pirho$Q[1] > 4
    pirho_Q2_great <- fitted_bisbmpop_pirho$Q[2] > 4

    Q1_great <- c(NA, NA, pi_Q1_great, rho_Q1_great, pirho_Q1_great)
    Q2_great <- c(NA, NA, pi_Q2_great, rho_Q2_great, pirho_Q2_great)

    ## Equals
    pi_Q1_equal <- fitted_bisbmpop_pi$Q[1] == 4
    pi_Q2_equal <- fitted_bisbmpop_pi$Q[2] == 4

    rho_Q1_equal <- fitted_bisbmpop_rho$Q[1] == 4
    rho_Q2_equal <- fitted_bisbmpop_rho$Q[2] == 4

    pirho_Q1_equal <- fitted_bisbmpop_pirho$Q[1] == 4
    pirho_Q2_equal <- fitted_bisbmpop_pirho$Q[2] == 4

    Q1_equal <- c(NA, NA, pi_Q1_equal, rho_Q1_equal, pirho_Q1_equal)
    Q2_equal <- c(NA, NA, pi_Q2_equal, rho_Q2_equal, pirho_Q2_equal)

    data_frame_output <- data.frame(
        epsilon_alpha = rep(ea, 5),
        pi1 = matrix(pi1, nrow = 1),
        rho2 = matrix(rho2, nrow = 1),
        model = model,
        BICL = BICLs,
        Q1_less = Q1_less,
        Q1_equal = Q1_equal,
        Q1_great = Q1_great,
        Q2_less = Q2_less,
        Q2_equal = Q2_equal,
        Q2_great = Q2_great,
        row_ARI = row_ARIs,
        col_ARI = col_ARIs,
        repetition = rep(conditions[c,4], 5)
    )

    return(data_frame_output)
},
mc.cores = parallel::detectCores() - 1)

saveRDS(results, "./simulation/data")
devtools::load_all()
require("sbm")
require("pROC")

set.seed(1234)

eps <- 0.05

M <- 3

# Defining parameters
nr <- 100
nc <- 150
pir <- c(0.5, 0.3, 0.2)
pic <- c(0.5, 0.3, 0.2)
alpha <- matrix(c(
  0.6, eps, eps,
  eps, 0.4, eps,
  eps, eps, 0.7
), byrow = TRUE, nrow = length(pir), ncol = length(pic))

# Collections
collections <- list(
  iid = generate_bipartite_collection(nr, nc,
    pir, pic,
    alpha, M,
    model = "iid",
    return_memberships = TRUE
  ),
  pi = generate_bipartite_collection(nr, nc,
    pir, pic,
    alpha, M,
    model = "pi",
    return_memberships = TRUE
  ),
  rho = generate_bipartite_collection(nr, nc,
    pir, pic,
    alpha, M,
    model = "rho",
    return_memberships = TRUE
  ),
  pirho = generate_bipartite_collection(nr, nc,
    pir, pic,
    alpha, M,
    model = "pirho",
    return_memberships = TRUE
  )
)


conditions <- expand.grid(
  prop_NAs = seq(from = 0, to = 0.9, by = 0.1),
  model = c("iid", "pi", "rho", "pirho"),
  repetition = seq.int(3)
)

result_dataframe <- do.call("rbind", bettermc::mclapply(seq_len(nrow(conditions)), function(current) {
  # Looping over conditions
  prop_NAs <- conditions[current, ]$prop_NAs
  model <- as.character(conditions[current, ]$model)
  bipartite_collection <- collections[[model]]

  # This is a list of the M incidence matrices
  bipartite_collection_incidence <- lapply(seq.int(M), function(m) {
    bipartite_collection[[m]]$incidence_matrix
  })

  # Sampling values to replace by NAs
  NAs_index <- sample(
    seq_len(length(bipartite_collection_incidence[[1]])),
    floor(prop_NAs * length(bipartite_collection_incidence[[1]]))
  )

  real_val_NAs <- bipartite_collection_incidence[[1]][NAs_index]
  bipartite_collection_incidence[[1]][NAs_index] <- NA
  NAs_coordinates <- which(is.na(bipartite_collection_incidence[[1]]),
    arr.ind = TRUE
  )

  Z <- lapply(seq.int(M), function(m) {
    list(
      bipartite_collection[[m]]$row_blockmemberships,
      bipartite_collection[[m]]$col_blockmemberships
    )
  })

  start_time <- Sys.time()
  mybisbmpop <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence, colsbm_model = model,
    global_opts = list(
      nb_cores = parallel::detectCores() - 1, verbosity = 0
    ), silent_parallelization = TRUE
  )
  stop_time <- Sys.time()

  baseline_LBM <- estimate_colBiSBM(
    netlist = bipartite_collection_incidence[[1]], colsbm_model = model,
    global_opts = list(
      nb_cores = parallel::detectCores() - 1, verbosity = 0
    ), silent_parallelization = TRUE
  )

  # Predicted links
  X_hat <- mybisbmpop$best_fit$tau[[1]][[1]] %*% mybisbmpop$best_fit$alpha %*% t(mybisbmpop$best_fit$tau[[1]][[2]])

  # Compute ROC and AUC
  auc <- auc(real_val_NAs, X_hat[NAs_index])

  # Computing ARI on the NAs
  return(data.frame(
    prop_NAs = prop_NAs,
    model = model,
    repetition = conditions[current, ]$repetition,
    auc = auc,
    LBM_ari_row = aricode::ARI(
      Z[[1]][[1]],
      baseline_LBM$best_fit$Z[[1]][[1]]
    ),
    LBM_ari_col = aricode::ARI(
      Z[[1]][[2]],
      baseline_LBM$best_fit$Z[[1]][[2]]
    ),
    NA_net_ari_row = aricode::ARI(
      Z[[1]][[1]],
      mybisbmpop$best_fit$Z[[1]][[1]]
    ),
    NA_net_ari_col = aricode::ARI(
      Z[[1]][[2]],
      mybisbmpop$best_fit$Z[[1]][[2]]
    ),
    elapsed_secs = difftime(stop_time, start_time, units = "sec")
  ))
},
mc.cores = parallel::detectCores() - 1,
mc.progress = TRUE
))

saveRDS(
  result_dataframe,
  paste0(
    "simulation/data/",
    "NA_robustness_results-", format(Sys.time(), "%d-%m-%y_%H-%M"), ".Rds"
  )
)

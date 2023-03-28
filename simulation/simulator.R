# Source necessary files to work
## Parallelization
require("bettermc")
require("ggplot2")
require("ggnewscale")
source("R/utils.R")
source("R/R6class-fitBipartiteSBMPop.R")
source("R/R6class-lbmpop.R")

combine_matrices_print <- function(matrix_print1, matrix_print2, sep = "") {
    list1 <- strsplit(matrix_print1, "\n")[[1]]
    list2 <- strsplit(matrix_print2, "\n")[[1]]

    sep <- append(c(ifelse(is.null(sep), "", paste(rep(" ", length(list1) + nchar(sep) - 4),
        collapse = ""
    ))), rep(sep, length(list1) - 1))

    return(paste(as.list(paste(list1, sep, list2)), collapse = "\n"))
}

eps <- 0.05

nr <- 250
nc <- 250

pir <- as.vector(gtools::rdirichlet(1, c(8, 5, 2)))
pic <- as.vector(gtools::rdirichlet(1, c(8, 5, 2)))

Q <- c(length(pir), length(pic))

isParallelized <- TRUE
first_alpha <- matrix(
    c(
        0.9, eps, eps,
        eps, 0.8, eps,
        eps, eps, 0.9
    ),
    byrow = TRUE,
    nrow = Q[1],
    ncol = Q[2]
)

second_alpha <- matrix(
    c(
        0.9, 0.8, 0.4,
        0.8, 0.2, eps,
        0.4, eps, eps
    ),
    byrow = TRUE,
    nrow = Q[1],
    ncol = Q[2]
)

filename <- "divergence_modular_to_nested"
filename <- paste0(getwd(), "/simulation/data/", filename, "_", format(Sys.time(), "%d-%m-%y_%X"), ".Rds")


complete_tibble <- NULL

# Blur goes from 0 to 0.9 by 0.1 step
# M takes 1,2 and 5
condition_matrix <- expand.grid(blur = seq(0, 1, by = 0.1), M = c(1, 2, 5))

if (isParallelized) {
    n.cores <- parallel::detectCores() - 1
} else {
    n.cores <- 1L
}

results <- bettermc::mclapply(seq.int(nrow(condition_matrix)), function(condition_row) {
    divergence_parameter <- condition_matrix[condition_row, 1]
    M <- condition_matrix[condition_row, 2]

    first_structure_collection <- generate_bipartite_collection(nr, nc, pir, pic, first_alpha, M)


    # Diverging alpha
    diverging_alpha <- (1 - divergence_parameter) * first_alpha + divergence_parameter * second_alpha
    # The second structure will receive a slowly diverging alpha
    diverging_collection <- generate_bipartite_collection(nr, nc, pir, pic, diverging_alpha, M)

    bipartite_collection <- append(first_structure_collection, diverging_collection)

    collection_incidence <- lapply(seq_along(bipartite_collection), function(m) {
        bipartite_collection[[m]]$incidence_matrix
    })

    collection_clustering <- lapply(seq_along(bipartite_collection), function(m) {
        list(
            bipartite_collection[[m]]$row_clustering,
            bipartite_collection[[m]]$col_clustering
        )
    })

    current_lbmpop <- lbmpop$new(
        netlist = collection_incidence,
        model = "bernoulli",
        free_mixture = FALSE,
        free_density = FALSE,
        global_opts = list(
            verbosity = 0,
            plot_details = 0,
            nb_cores = 1
        )
    )

    current_lbmpop$optimize()

    current_tibble <- dplyr::bind_rows(lapply(seq_along(current_lbmpop$best_fit$MAP$Z), function(m) {
        tibble::tibble(
            Current_M = M,
            Network_id = current_lbmpop$best_fit$net_id[[m]],
            row_ARI = aricode::ARI(
                as.vector(current_lbmpop$best_fit$MAP$Z[[m]][[1]]),
                collection_clustering[[m]][[1]]
            ),
            col_ARI = aricode::ARI(
                as.vector(current_lbmpop$best_fit$MAP$Z[[m]][[2]]),
                collection_clustering[[m]][[2]]
            ),
            divergence = divergence_parameter,
            sep_LBM_BICL = sum(current_lbmpop$sep_LBM_BICL),
            BICL = current_lbmpop$best_fit$BICL
        )
    }))
    current_tibble
}, mc.cores = n.cores)

complete_tibble <- dplyr::bind_rows(results)

data_to_save <- list(
    pir = pir,
    pic = pic,
    first_alpha = first_alpha,
    second_alpha = second_alpha,
    results = complete_tibble
)

saveRDS(data_to_save, file = filename)

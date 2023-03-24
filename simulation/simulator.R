# Source necessary files to work
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

real_alpha <- matrix(
    c(
        0.9, eps, eps,
        eps, 0.8, eps,
        eps, eps, 0.9
    ),
    byrow = TRUE,
    nrow = Q[1],
    ncol = Q[2]
)

filename <- "modular"
filename <- paste0(getwd(), "/simulation/data/", filename, "_", format(Sys.time(), "%d-%m-%y_%X"), ".Rds")


# Function
real_alpha_print <- paste0("Real alpha",
    capture.output(print(real_alpha)),
    collapse = "\n"
)

er_probability <- 0.25
er_alpha <- matrix(rep(er_probability, Q[1] * Q[2]), nrow = Q[1], ncol = Q[2])

er_alpha_print <- paste0(capture.output(print(er_alpha)), collapse = "\n")
complete_tibble <- NULL

# Blur goes from 0 to 0.9 by 0.1 step
# For M from 2 to 10 by one increase
condition_matrix <- expand.grid(seq(0, 1, by = 0.1), seq(2, 10, by = 1))
# Progressbar
progress_bar <- progress::progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = nrow(condition_matrix),
    complete = "=", # Completion bar character
    incomplete = "-", # Incomplete bar character
    current = ">", # Current bar character
    clear = FALSE, # If TRUE, clears the bar when finish
    width = 100, # Width of the progress bar
    show_after = 0.0
)
for (condition_row in 1:nrow(condition_matrix)) {

    blur_parameter <- condition_matrix[condition_row, 1]
    M <- condition_matrix[condition_row, 2]

    alpha <- (1 - blur_parameter) * real_alpha + blur_parameter * er_alpha
    alpha_print <- paste0(capture.output(print(alpha)), collapse = "\n")

    bipartite_collection <- generate_bipartite_collection(nr, nc, pir, pic, alpha, M)
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
            verbosity = 4,
            plot_details = 0
        )
    )

    current_lbmpop$optimize()

    current_tibble <- dplyr::bind_rows(lapply(seq_along(current_lbmpop$best_fit$MAP$Z), function(m) {
        tibble::tibble(
            Current_M = M,
            Network_id = current_lbmpop$best_fit$net_id[[m]],
            row_ARI = aricode::ARI(
                current_lbmpop$best_fit$MAP$Z[[m]][[1]],
                collection_clustering[[m]][[1]]
            ),
            col_ARI = aricode::ARI(
                current_lbmpop$best_fit$MAP$Z[[m]][[2]],
                collection_clustering[[m]][[2]]
            ),
            blur = blur_parameter
        )
    }))
    complete_tibble <- dplyr::bind_rows(complete_tibble, current_tibble)
    progress_bar$tick()

}

data_to_save <- list(
    pir = pir,
    pic = pic,
    real_alpha = real_alpha,
    er_probability = er_probability,
    results = complete_tibble
)

saveRDS(data_to_save, file = filename)

# fitSimpleSBMPop$set(
#   "public", "plot",

# )

# bmpop$set(
#   "public",
#   "plot",

# )

# ## Bipartite plots
# fitBipartiteSBMPop$set(
#   "public", "plot",

# )


#' Plot matrix summaries of the collection mesoscale structure
#'
#' @param x a fitSimpleSBMPOP object.
#' @param type The type of the plot. Could be "graphon", "meso" or "block".
#' @param ord A reordering of the blocks.
#' @param mixture Should the block proportions of each network be plotted as
#' well?
#' @param net_id Use to plot only on network in "graphon" view.
#' @param ... Further argument to be passed
#' @return A plot, a ggplot2 object.
#' @export
#'
#' @examples
#' # Trivial example with Gnp networks:
#' Net <- lapply(
#'   list(.7, .7, .2, .2),
#'   function(p) {
#'     A <- matrix(0, 15, 15)
#'     A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1
#'     A <- A + t(A)
#'   }
#' )
#' \dontrun{
#' cl <- estimate_colSBM(Net,
#'   colsbm_model = "delta",
#'   directed = FALSE,
#'   distribution = "bernoulli",
#'   nb_run = 1
#' )
#' plot(cl$best_fit)
#' }
plot.fitSimpleSBMPop <- function(x, type = "graphon",
                                 ord = NULL, mixture = FALSE, net_id = 1, ...) {
  stopifnot(inherits(x, "fitSimpleSBMPop"))
  p <- x$plot(
    type = type, ord = ord, mixture = mixture,
    net_id = net_id, ...
  )
  p
}

# plot_colsbm_nomix <- function(fit, ord = NULL,  title = NULL, tag = NULL) {
#   if (fit$free_mixture) {pim <- fit$pi} else {pim <- fit$pim}
#   #ord <- order(colSums(fit$alpha > .01) - rowSums(fit$alpha >.01)/2, decreasing = FALSE)
#   if(is.null(ord)) ord <- order(diag(fit$alpha), decreasing = TRUE)
#   p_alpha_pr_res <- fit$alpha[ord, ord] %>% t() %>%
#     reshape2::melt() %>%
#     ggplot(aes(x = Var1, y = Var2, fill = value)) +
#     geom_tile() +
#     scale_fill_gradient2("alpha", low = "white", high = "red") +
#     geom_hline(yintercept =  seq(fit$Q)+.5) +
#     geom_vline(xintercept =  seq(fit$Q)+.5) +
#     scale_y_reverse() +
#     theme_bw(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
#     xlab("") + ylab("") + coord_fixed(expand = FALSE) +
#     labs(tag =  tag)
#   #  scale_y_reverse()
#   if (fit$free_density) {
#     xl <- paste(round(fit$delta,1))
#   } else {
#     xl <-  ""
#   }
#   df <-     map_dfc(seq_along(fit$net_id) ,
#                     function(m)  fit$pim[[m]][ord])
#   names(df) <- fit$net_id
#   p_pi_pi_pr_res <-
#     df %>%
#     #    rename() %>%
#     mutate(q = seq(fit$Q)) %>%
#     pivot_longer(cols = -c(q)) %>%
#     ggplot(aes(fill = as.factor(q), y = name, x = value)) +
#     geom_col() +
#     coord_flip(expand = FALSE) +
#     scale_fill_brewer("Block", type = "qual", palette = "Paired", direction = -1) +
#     ylab("") +
#     ylab(xl) +
#     theme_bw(base_size = 15)
#
#   ( p_alpha_pr_res + p_pi_pi_pr_res ) + patchwork::plot_layout(guides = 'collect', widths = c(.65,.35)) +
#     plot_annotation(title = NULL)
# }
# A_priest[order(sbm_priest$best_fit$Z[[1]]),
#          order(sbm_priest$best_fit$Z[[1]])] %>%
#   reshape2::melt() %>%
#   ggplot(aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile(show.legend = FALSE) +
#   geom_hline(yintercept = cumsum(tabulate(sbm_priest$best_fit$Z[[1]])[1:5])+.5,
#              col = "red", size = .5) +
#   geom_vline(xintercept = cumsum(tabulate(sbm_priest$best_fit$Z[[1]])[1:5])+.5,
#              col = "red", size = .5) +
#   scale_fill_gradient(low = "white", high = "black") +
#   ylab("") + xlab("Priests") +
#   scale_x_discrete(#limits = rev,
#     breaks = "") +
#   scale_y_reverse() +
#   # scale_y_discrete(limits = rev, breaks = "",#label = rev(custom_lab3),
#   #                  guide = guide_axis(angle = 0) ) +
#   coord_equal(expand = FALSE) +
#   theme_bw(base_size = 15) +
#   theme(axis.ticks =  element_blank())

#' Plot the trace of the different criteria in function of the number of
#' clusters
#'
#' @param x a bmpop object.
#' @param type The type of the plot. "trace".
#' @param ... Further argument to be passed
#' @return A plot, a ggplot2 object.
#' @export
#'
#'
plot.bmpop <- function(x, type = "trace", ...) {
  stopifnot(inherits(x, "bmpop"))
  p <- x$plot(type = type, ...)
  p
}

#' Plot matrix summaries of the collection mesoscale structure
#'
#' @param x a fitBipartiteSBMPop object.
#' @param type The type of the plot. Could be "graphon", "meso" or "block".
#' @param oRow A reordering of the row blocks.
#' @param oCol A reordering of the column blocks.
#' @param mixture Should the block proportions of each network be plotted as
#' well?
#' @param net_id Use to plot only on network in "graphon" view.
#' @param ... Further argument to be passed
#' @return A plot, a ggplot2 object.
#' @export
#'
#' @examples
#' alpha1 <- matrix(c(0.8, 0.1, 0.2, 0.7), byrow = TRUE, nrow = 2)
#'
#' first_collection <- generate_bipartite_collection(
#'   nr = 50, nc = 25,
#'   pi = c(0.5, 0.5), rho = c(0.5, 0.5),
#'   alpha = alpha1, M = 2
#' )
#'
#' \dontrun{
#' # A collection where joint modelisation makes sense
#' cl_joint <- estimate_colBiSBM(
#'   netlist = first_collection,
#'   colsbm_model = "iid",
#'   global_opts = list(nb_cores = parallelly::availableCores(omit = 1L))
#' )
#' plot(cl_joint$best_fit)
#' }
plot.fitBipartiteSBMPop <- function(x, type = "graphon", oRow = NULL, oCol = NULL, mixture = FALSE, net_id = 1, ...) {
  stopifnot(inherits(x, "fitBipartiteSBMPop"))
  p <- x$plot(
    type = type, oRow = oRow, oCol = oCol, mixture = mixture,
    net_id = net_id, ...
  )
  p
}

#' Plot the state-space exploration plot for a bipartite collection object
#'
#' @param x a bisbmpop object.
#' @param ... other arguments to pass to the plot.
#'
#' @return A plot, a ggplot2 object.
#' @export
plot.bisbmpop <- function(x, ...) {
  stopifnot(inherits(x, "bisbmpop"))
  p <- x$plot(...)
  p
}
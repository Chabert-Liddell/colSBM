fitSimpleSBMPop$set("public", "plot",
#' Title
#'
#' @importFrom patchwork
#' @importFrom reshape2
#' @importFrom purrr
  function(type = "graphon", ord = NULL, mixture = FALSE, net_id = NULL, ...) {
  if(is.null(ord)) ord <- order(diag(self$alpha), decreasing = TRUE)
  p <- switch(
    type,
    graphon = {
      (self$alpha[ord, ord]*mean(self$delta)) %>% t() %>%
        reshape2::melt() %>%
        dplyr::mutate(xmax = rep(c(0,cumsum(self$pi[[net_id]][ord][1:(self$Q-1)])), self$Q),
                      xmin = rep(cumsum(self$pi[[net_id]][ord]), self$Q),
                      ymax = rep(c(0,cumsum(self$pi[[net_id]][ord][1:(self$Q-1)])), each = self$Q),
                      ymin = rep(cumsum(self$pi[[net_id]][ord]), each = self$Q)) %>%
        ggplot2::ggplot(ggplot2::aes(xmin = xmin, ymin = ymin,
                                     xmax = xmax, ymax = ymax, fill = value)) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", mid = "red", midpoint = 1) +
        ggplot2::geom_hline(yintercept = cumsum(self$pi[[net_id]][ord][1:(self$Q-1)]), size = .2) +
        ggplot2::geom_vline(xintercept = cumsum(self$pi[[net_id]][ord][1:(self$Q-1)]), size = .2) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_equal(expand = FALSE)
    },
    meso = {
      if (self$free_mixture) {pim <- self$pi} else {pim <- self$pim}
      p_alpha <- self$alpha[ord, ord] %>% t() %>%
        reshape2::melt() %>%
        ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", high = "red") +
        ggplot2::geom_hline(yintercept =  seq(self$Q)+.5) +
        ggplot2::geom_vline(xintercept =  seq(self$Q)+.5) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
        ggplot2::xlab("") + ggplot2::ylab("") +
        ggplot2::coord_fixed(expand = FALSE)
      #  scale_y_reverse()
      if (self$free_density) {
        xl <- paste(round(self$delta,1))
      } else {
        xl <-  ""
      }
      # df <-     purrr::map_dfc(seq_along(self$net_id) ,
      #                   function(m)  self$pim[[m]][ord])
      df <-  purrr::list_cbind(
        lapply(seq_along(self$net_id),
               function(m)  tibble::tibble(self$pim[[m]][ord])),
        name_repair =  "universal_quiet")
      names(df) <- self$net_id
      if (mixture) {
        p_pi <-
          df %>%
          #    rename() %>%
          dplyr::mutate(q = seq(self$Q)) %>%
          tidyr::pivot_longer(cols = -c(q)) %>%
          ggplot2::ggplot(ggplot2::aes(fill = as.factor(q), y = name, x = value)) +
          ggplot2::geom_col() +
          ggplot2::coord_flip(expand = FALSE) +
          ggplot2::scale_fill_brewer("Block", type = "qual", palette = "Paired",
                                     direction = -1) +
          ggplot2::ylab("") +
          ggplot2::xlab(xl) +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 70,
                                                             vjust = 0.5))
        p_alpha <- ( p_alpha | p_pi) +
          patchwork::plot_layout(guides = 'collect', widths = c(.65,.35)) +
          patchwork::plot_annotation(title = NULL)
      }
      return(p_alpha)
    },
    "block" = {as.matrix(self$A[[net_id]])[order(self$Z[[net_id]]),
                               order(self$Z[[net_id]])] %>%
      reshape2::melt() %>%
      dplyr::mutate(con = self$alpha[self$Z[[net_id]], self$Z[[net_id]]][
        order(ord[self$Z[[net_id]]]),
        order(ord[self$Z[[net_id]]])] %>%
                 reshape2::melt() %>%
                 dplyr::pull(value)) %>%
      ggplot2::ggplot(ggplot2::aes(x = Var2, y = Var1,
                                   fill = value, alpha = value)) +
      ggplot2::geom_tile(ggplot2::aes(alpha = con),
                         fill = "red", size = 0, show.legend = FALSE) +
      ggplot2::geom_tile(show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = cumsum(tabulate(self$Z[[net_id]])[1:(self$Q-1)])+.5,
                 col = "red", linewidth = .5) +
        ggplot2::geom_vline(xintercept = cumsum(tabulate(self$Z[[net_id]])[(self$Q):2])+.5,
                 col = "red", linewidth = .5) +
        ggplot2::scale_fill_gradient(low = "white", high = "black") +
        ggplot2::ylab("") + ggplot2::xlab(self$net_id[net_id]) +
        ggplot2::scale_x_discrete(limits = rev,
        breaks = "") +
        #ggplot2::scale_y_reverse() +
       ggplot2::scale_y_discrete(#limits = rev,
                                 breaks = "",#label = rev(custom_lab3),
                        guide = ggplot2::guide_axis(angle = 0) ) +
        ggplot2::scale_alpha(range = c(0,1)) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(axis.ticks =  ggplot2::element_blank())}
  )
  return(p)
})


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
plot.fitSimpleSBMPop <- function(x, type = "graphon",
                                 ord = NULL, mixture = FALSE, net_id = 1, ...) {
  stopifnot(inherits(x, "fitSimpleSBMPop"))
  p <- x$plot(type = type, ord = ord, mixture = mixture,
              net_id = net_id, ...)
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



bmpop$set(
  "public",
  "plot",
#' Title
#'
#' @param type
#' @param ...
#'
#' @importFrom tibble
#' @importFrom ggplot2
#' @importFrom tidyr
#' @importFrom dplyr
  function(type = "trace", ...) {
    tb <-           tibble::tibble(
      Q = seq(length(self$BICL)),
      ICL = self$ICL,
      BICL = self$BICL,
      vbound = self$vbound)
  if (! is.null(self$ICL_sbm)) {
    tb %>% dplyr::mutate(SBM = self$ICL_sbm)
  }
    tb %>%
      tidyr::pivot_longer(cols = -Q, names_to = "Criterion") %>%
      ggplot2::ggplot(ggplot2::aes(x = Q, y = value,
                                   linetype = Criterion, color = Criterion,
                                   shape = Criterion)) +
      ggplot2::annotate(geom = "rect",
                        xmin = max(which.max(self$BICL) -
                                     self$global_opts$depth, 1),
                        xmax = min(which.max(self$BICL) +
                                     self$global_opts$depth, length(self$BICL)),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "gray90"
      ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 3) +
      ggplot2::ylab("") +
      ggplot2::theme_bw()
  }
)

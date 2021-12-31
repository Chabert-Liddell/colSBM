#' Title
#'
#' @param fit
#' @param ord
#' @param title
#' @param tag
#'
#' @return
#' @export
#'
#' @examples
plot.fitSimpleSBMPop <- function(x, ord = NULL,  title = NULL, tag = NULL) {
  stopifnot(inherits(x, "fitSimpleSBMPop"))
  p <- x$plot(ord = ord, title = title, tag = tag)
  p
}

plot_colsbm_nomix <- function(fit, ord = NULL,  title = NULL, tag = NULL) {
  if (fit$free_mixture) {pim <- fit$pi} else {pim <- fit$pim}
  #ord <- order(colSums(fit$alpha > .01) - rowSums(fit$alpha >.01)/2, decreasing = FALSE)
  if(is.null(ord)) ord <- order(diag(fit$alpha), decreasing = TRUE)
  p_alpha_pr_res <- fit$alpha[ord, ord] %>% t() %>%
    reshape2::melt() %>%
    ggplot(aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2("alpha", low = "white", high = "red") +
    geom_hline(yintercept =  seq(fit$Q)+.5) +
    geom_vline(xintercept =  seq(fit$Q)+.5) +
    scale_y_reverse() +
    theme_bw(base_size = 15, base_rect_size = 1, base_line_size  = 1) +
    xlab("") + ylab("") + coord_fixed(expand = FALSE) +
    labs(tag =  tag)
  #  scale_y_reverse()
  if (fit$free_density) {
    xl <- paste(round(fit$delta,1))
  } else {
    xl <-  ""
  }
  df <-     map_dfc(seq_along(fit$net_id) ,
                    function(m)  fit$pim[[m]][ord])
  names(df) <- fit$net_id
  p_pi_pi_pr_res <-
    df %>%
    #    rename() %>%
    mutate(q = seq(fit$Q)) %>%
    pivot_longer(cols = -c(q)) %>%
    ggplot(aes(fill = as.factor(q), y = name, x = value)) +
    geom_col() +
    coord_flip(expand = FALSE) +
    scale_fill_brewer("Block", type = "qual", palette = "Paired", direction = -1) +
    ylab("") +
    ylab(xl) +
    theme_bw(base_size = 15)

  ( p_alpha_pr_res + p_pi_pi_pr_res ) + patchwork::plot_layout(guides = 'collect', widths = c(.65,.35)) +
    plot_annotation(title = NULL)
}
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

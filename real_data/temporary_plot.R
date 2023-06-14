#' The function to plot the fitBipartite objects
#' @importFrom patchwork
#' @importFrom reshape2
#' @importFrom purrr
plotFitBipartite <- function(model, type = "graphon", oRow = NULL, oCol = NULL, mixture = FALSE, net_id = NULL, ...) {
  # The below order use mean over all networks to have a consistent display
  if (is.null(oRow)) {
    if (model$Q[2] == 1) {
      mean_rho <- 1
    } else {
      mean_rho <- matrixStats::rowMeans2(sapply(model$pi, function(pi) pi[[2]]))
    }
    oRow <- order(model$alpha %*% mean_rho, decreasing = TRUE)
  }
  if (is.null(oCol)) {
    if (model$Q[1] == 1) {
      mean_pi <- 1
    } else {
      mean_pi <- matrixStats::rowMeans2(sapply(model$pi, function(pi) pi[[1]]))
    }
    oCol <- order(mean_pi %*% model$alpha, decreasing = TRUE)
  }
  p <- switch(type,
    graphon = {
      (model$alpha[oRow, oCol] * mean(model$delta)) %>%
        t() %>%
        reshape2::melt() %>%
        dplyr::mutate(
          xmin = rep(c(0, cumsum(model$pi[[net_id]][[2]][oCol][1:(model$Q[2] - 1)])), model$Q[1]),
          ymin = rep(c(0, cumsum(model$pi[[net_id]][[1]][oRow][1:(model$Q[1] - 1)])), each = model$Q[2]),
          xmax = rep(cumsum(model$pi[[net_id]][[2]][oCol]), model$Q[1]),
          ymax = rep(cumsum(model$pi[[net_id]][[1]][oRow]), each = model$Q[2])
        ) %>%
        ggplot2::ggplot(ggplot2::aes(
          xmin = xmin, ymin = ymin,
          xmax = xmax, ymax = ymax, fill = value
        )) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", mid = "red", midpoint = 1) +
        ggplot2::geom_hline(yintercept = cumsum(model$pi[[net_id]][[1]][oRow][1:(model$Q[1] - 1)]), size = .2) +
        ggplot2::geom_vline(xintercept = cumsum(model$pi[[net_id]][[2]][oCol][1:(model$Q[2] - 1)]), size = .2) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("Column Blocks") +
        ggplot2::ylab("Row Blocks") +
        ggplot2::coord_equal(expand = FALSE)
    },
    meso = {
      p_alpha <- model$alpha[oRow, oCol] %>%
        t() %>%
        reshape2::melt() %>%
        ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", high = "red") +
        ggplot2::geom_hline(yintercept = seq(model$Q[1]) + .5) +
        ggplot2::geom_vline(xintercept = seq(model$Q[2]) + .5) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_fixed(expand = FALSE)
      #  scale_y_reverse()
      if (model$free_density) {
        xl <- paste(round(model$delta, 1))
      } else {
        xl <- ""
      }
      df_pi <- purrr::map_dfc(
        seq_along(model$net_id),
        function(m) setNames(data.frame(model$pim[[m]][[1]][oRow]), m)
      )
      df_rho <- purrr::map_dfc(
        seq_along(model$net_id),
        function(m) setNames(data.frame(model$pim[[m]][[2]][oCol]), m)
      )
      # names(df_pi) <- model$net_id
      if (mixture) {
        p_pi <-
          df_pi %>%
          #    rename() %>%
          dplyr::mutate(q = seq(model$Q[1])) %>%
          tidyr::pivot_longer(cols = -c(q)) %>%
          mutate(Proportion = value) %>%
          ggplot2::ggplot(ggplot2::aes(
            fill = as.factor(q), y = name,
            x = Proportion
          )) +
          ggplot2::geom_col() +
          ggplot2::coord_flip(expand = FALSE) +
          ggplot2::scale_fill_brewer("Row block",
            type = "qual", palette = "Paired",
            direction = -1
          ) +
          ggplot2::guides(fill = ggplot2::guide_legend(
            ncol = model$Q[1] %/% 3 + 1,
            byrow = TRUE
          )) +
          ggplot2::ylab("") +
          ggplot2::ylab(xl) +
          ggplot2::theme_bw(base_size = 15)
        p_rho <- df_rho %>%
          #    rename() %>%
          dplyr::mutate(q = seq(model$Q[2])) %>%
          tidyr::pivot_longer(cols = -c(q)) %>%
          mutate(Proportion = value) %>%
          ggplot2::ggplot(ggplot2::aes(
            fill = as.factor(q), y = name,
            x = Proportion
          )) +
          ggplot2::geom_col() +
          # ggplot2::coord_flip(expand = FALSE) +
          ggplot2::scale_fill_brewer("Column block",
            type = "qual", palette = "Set2",
            direction = -1
          ) +
          ggplot2::guides(fill = ggplot2::guide_legend(
            ncol = model$Q[2] %/% 3 + 1,
            byrow = TRUE
          )) +
          ggplot2::ylab("") +
          ggplot2::ylab(xl) +
          ggplot2::theme_bw(base_size = 15)
        # Merging the plots with patchwork
        mixture_layout <- "
                            ##CCCC
                            ##CCCC
                            RRAAAA
                            RRAAAA
                            RRAAAA
                            "
        p_alpha <- patchwork::wrap_plots(
          R = p_pi, C = p_rho, A = p_alpha,
          design = mixture_layout
        ) +
          patchwork::plot_layout(
            guides = "collect",
            design = mixture_layout
          )
      }
      return(p_alpha)
    },
    "block" = {
      as.matrix(model$A[[net_id]])[
        order(model$Z[[net_id]][[1]]),
        order(model$Z[[net_id]][[2]])
      ] %>%
        reshape2::melt() %>%
        ggplot2::ggplot(ggplot2::aes(x = Var2, y = rev(Var1), fill = value)) +
        ggplot2::geom_tile(show.legend = FALSE) +
        ggplot2::geom_hline(
          yintercept = cumsum(tabulate(model$Z[[net_id]][[1]])[model$Q[1]:2]) + .5,
          col = "red", size = .5
        ) +
        ggplot2::geom_vline(
          xintercept = cumsum(tabulate(model$Z[[net_id]][[2]])[1:(model$Q[2] - 1)]) + .5,
          col = "red", size = .5
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = "black", na.value = "transparent") +
        ggplot2::ylab("") +
        ggplot2::xlab(model$net_id[net_id]) +
        ggplot2::scale_x_discrete(
          breaks = ""
        ) +
        # ggplot2::scale_y_reverse() +
        ggplot2::scale_y_discrete(
          breaks = "",
          guide = ggplot2::guide_axis(angle = 0)
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank())
    }
  )
  return(p)
}

#' Plot matrix summaries of the collection mesoscale structure
#'
#' @param x a fitBipartiteSBMPop object.
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
plot.fitBipartiteSBMPop <- function(x, type = "graphon", oRow = NULL, oCol = NULL, mixture = FALSE, net_id = 1, ...) {
  stopifnot(inherits(x, "fitBipartiteSBMPop"))
  p <- plotFitBipartite(x,
    type = type, oRow = oRow, oCol = oCol, mixture = mixture,
    net_id = net_id, ...
  )
  p
}

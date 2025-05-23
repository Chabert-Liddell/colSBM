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
plot.fitSimpleSBMPop <- function(
    x,
    type = "graphon",
    ord = NULL,
    mixture = FALSE,
    net_id = 1L,
    ...) {
  if (is.null(ord)) ord <- order(diag(x$alpha), decreasing = TRUE)
  p <- switch(type,
    graphon = {
      if (x$Q == 1) {
        ymin <- rep(0, each = x$Q)
        ymax <- rep(1, each = x$Q)
        xmin <- rep(0, x$Q)
        xmax <- rep(1, x$Q)
      } else {
        xmax <- rep(c(0, cumsum(x$pi[[net_id]][ord][1:(x$Q - 1)])), x$Q)
        xmin <- rep(cumsum(x$pi[[net_id]][ord]), x$Q)
        ymax <- rep(c(0, cumsum(x$pi[[net_id]][ord][1:(x$Q - 1)])), each = x$Q)
        ymin <- rep(cumsum(x$pi[[net_id]][ord]), each = x$Q)
      }
      (x$alpha[ord, ord] * mean(x$delta)) %>%
        t() %>%
        reshape2::melt() %>%
        dplyr::mutate(
          xmax = xmax,
          xmin = xmin,
          ymax = ymax,
          ymin = ymin
        ) %>%
        ggplot2::ggplot(ggplot2::aes(
          xmin = xmin, ymin = ymin,
          xmax = xmax, ymax = ymax, fill = value
        )) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", mid = "red", midpoint = 1) +
        ggplot2::geom_hline(yintercept = cumsum(x$pi[[net_id]][ord][1:(x$Q - 1)]), linewidth = .2) +
        ggplot2::geom_vline(xintercept = cumsum(x$pi[[net_id]][ord][1:(x$Q - 1)]), linewidth = .2) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_equal(expand = FALSE)
    },
    meso = {
      if (x$free_mixture) {
        pim <- x$pi
      } else {
        pim <- x$pim
      }
      p_alpha <- x$alpha[ord, ord] %>%
        t() %>%
        reshape2::melt() %>%
        ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2("alpha", low = "white", high = "red") +
        ggplot2::geom_hline(yintercept = seq(x$Q) + .5) +
        ggplot2::geom_vline(xintercept = seq(x$Q) + .5) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_fixed(expand = FALSE)
      #  scale_y_reverse()
      if (x$free_density) {
        xl <- paste(round(x$delta, 1))
      } else {
        xl <- ""
      }
      df <- purrr::map_dfc(
        seq_along(x$net_id),
        function(m) {
          setNames(
            data.frame(x$pim[[m]][ord]),
            x$net_id[m]
          )
        }
      )
      names(df) <- x$net_id
      if (mixture) {
        p_pi <-
          df %>%
          #    rename() %>%
          dplyr::mutate(q = seq(x$Q)) %>%
          tidyr::pivot_longer(cols = -c(q)) %>%
          ggplot2::ggplot(ggplot2::aes(fill = as.factor(q), y = name, x = value)) +
          ggplot2::geom_col() +
          ggplot2::coord_flip(expand = FALSE) +
          ggplot2::scale_fill_brewer("Block",
            type = "qual", palette = "Paired",
            direction = -1
          ) +
          ggplot2::ylab("") +
          ggplot2::ylab(xl) +
          ggplot2::theme_bw(base_size = 15)
        p_alpha <- (p_alpha | p_pi) + patchwork::plot_layout(guides = "collect", widths = c(.65, .35)) +
          patchwork::plot_annotation(title = NULL)
      }
      return(p_alpha)
    },
    "block" = {
      Z <- factor(x$Z[[net_id]], levels = rev(ord))
      as.matrix(x$A[[net_id]])[
        order(Z),
        order(Z)
      ] %>% # t() %>%
        reshape2::melt() %>%
        dplyr::mutate(con = x$alpha[x$Z[[net_id]], x$Z[[net_id]]][
          order(Z),
          order(Z)
        ] %>% # t() %>%#[x$Z[[net_id]], x$Z[[net_id]]][
          #  order(Z),
          #    order(Z)] %>%
          reshape2::melt() %>%
          dplyr::pull(value)) %>%
        ggplot2::ggplot(ggplot2::aes(x = as.factor(Var2), y = as.factor(Var1), fill = value, alpha = value)) +
        ggplot2::geom_tile(ggplot2::aes(alpha = con),
          fill = "red", linewidth = 0, show.legend = FALSE
        ) +
        ggplot2::geom_tile(show.legend = FALSE) +
        ggplot2::geom_hline(
          yintercept = cumsum(tabulate(Z)[1:(x$Q - 1)]) + .5,
          col = "red", linewidth = .5
        ) +
        ggplot2::geom_vline(
          xintercept = cumsum(tabulate(Z)[ifelse(x$Q > 1, (x$Q):2, 1)]) + .5,
          col = "red", linewidth = .5
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = "black") +
        ggplot2::ylab("") +
        ggplot2::xlab(x$net_id[net_id]) +
        ggplot2::scale_x_discrete(
          limits = rev,
          breaks = ""
        ) +
        # ggplot2::scale_y_reverse() +
        ggplot2::scale_y_discrete( # limits = rev,
          breaks = "", # label = rev(custom_lab3),
          guide = ggplot2::guide_axis(angle = 0)
        ) +
        ggplot2::scale_alpha(range = c(0, 1)) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank())
    }
  )
  return(p)
}

# plot_colsbm_nomix <- function(fit, ord = NULL,  title = NULL, tag = NULL) {
#   if (fit$free_mixture) {pim <- fit$pi} else {pim <- fit$pim}
#   #ord <- order(colSums(fit$alpha > .01) - rowSums(fit$alpha >.01)/2, decreasing = FALSE)
#   if(is.null(ord)) ord <- order(diag(fit$alpha), decreasing = TRUE)
#   p_alpha_pr_res <- fit$alpha[ord, ord] |> t() |>
#     reshape2::melt() |>
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
#     df |>
#     #    rename() |>
#     mutate(q = seq(fit$Q)) |>
#     pivot_longer(cols = -c(q)) |>
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
#          order(sbm_priest$best_fit$Z[[1]])] |>
#   reshape2::melt() |>
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
#'
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate pull
#'
#' @export
#'
#'
plot.bmpop <- function(x, type = "trace", ...) {
  stopifnot(inherits(x, "bmpop"))
  tb <- tibble::tibble(
    Q = seq_along(x$BICL),
    ICL = x$ICL,
    BICL = x$BICL,
    vbound = x$vbound
  )
  if (!is.null(x$ICL_sbm)) {
    tb %>% dplyr::mutate(SBM = x$ICL_sbm)
  }
  p <- tb %>%
    tidyr::pivot_longer(cols = -Q, names_to = "Criterion") %>%
    ggplot2::ggplot(ggplot2::aes(
      x = Q, y = value,
      linetype = Criterion, color = Criterion,
      shape = Criterion
    )) +
    ggplot2::annotate(
      geom = "rect",
      xmin = max(which.max(x$BICL) -
        x$global_opts$depth, 1),
      xmax = min(which.max(x$BICL) +
        x$global_opts$depth, length(x$BICL)),
      ymin = -Inf,
      ymax = Inf,
      fill = "gray90"
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::ylab("") +
    ggplot2::theme_bw()
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
#' @param values Wether or not to plot values on the alpha, pi and rho
#' representation.
#' @param values_min The minimum value to plot the value.
#' @param net_id Use to plot only on network in "graphon" view.
#' @param ... Further argument to be passed
#' @return A plot, a ggplot2 object.
#' @export
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom stats na.omit
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom reshape2 melt
#' @importFrom matrixStats rowMeans2
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
plot.fitBipartiteSBMPop <- function(
    x, type = "graphon",
    oRow = NULL,
    oCol = NULL,
    mixture = FALSE,
    net_id = 1L,
    values = FALSE,
    values_min = 0.1,
    ...) {
  stopifnot(inherits(x, "fitBipartiteSBMPop"))
  p <- switch(type,
    graphon = {
      # The below order use net_id parameters
      if (x$Q[2] == 1) {
        mean_rho <- 1
      } else {
        mean_rho <- x[["pi"]][[net_id]][[2]]
      }
      if (is.null(oRow)) {
        oRow <- order(x$alpham[[net_id]] %*% mean_rho, decreasing = TRUE)
      }
      #  Once order of tile in block will be fix, need to go back in nullity cond
      if (x$Q[1] == 1) {
        mean_pi <- 1
      } else {
        mean_pi <- x[["pi"]][[net_id]][[1]]
      }
      if (is.null(oCol)) {
        oCol <- order(mean_pi %*% x$alpham[[net_id]], decreasing = TRUE)
      }
      if (x$Q[1] == 1) {
        ymin <- rep(0, each = x$Q[2])
        ymax <- rep(1, each = x$Q[2])
      } else {
        ymin <- rep(c(0, cumsum(x$pi[[net_id]][[1]][oRow][1:(x$Q[1] - 1)])), each = x$Q[2])
        ymax <- rep(c(cumsum(x$pi[[net_id]][[1]][oRow])), each = x$Q[2])
      }
      if (x$Q[2] == 1) {
        xmin <- rep(0, x$Q[1])
        xmax <- rep(1, x$Q[1])
      } else {
        xmin <- rep(c(0, cumsum(x$pi[[net_id]][[2]][oCol][1:(x$Q[2] - 1)])), x$Q[1])
        xmax <- rep(cumsum(x$pi[[net_id]][[2]][oCol]), x$Q[1])
      }
      p_graphon <- (x$alpha[oRow, oCol]) |>
        t() |>
        reshape2::melt() |>
        dplyr::mutate(
          xmin = xmin,
          ymin = ymin,
          xmax = xmax,
          ymax = ymax
        ) |>
        ggplot2::ggplot(ggplot2::aes(
          xmin = xmin, ymin = ymin,
          xmax = xmax, ymax = ymax, fill = value
        )) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_gradient2("alpha",
          low = "white", mid = "red",
          midpoint = 1, limits = c(0, ifelse(x$distribution == "bernoulli", 1, max(x$alpha)))
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "α")) +
        ggplot2::geom_hline(yintercept = cumsum(x$pi[[net_id]][[1]][oRow][1:(x$Q[1] - 1)]), linewidth = .2) +
        ggplot2::geom_vline(xintercept = cumsum(x$pi[[net_id]][[2]][oCol][1:(x$Q[2] - 1)]), linewidth = .2)
      if (values) {
        p_graphon <- p_graphon +
          ggplot2::geom_text(ggplot2::aes(x = (Var2 - min(Var2)) / max(Var2), y = (Var1 - 0.5 * min(Var1)) / max(Var1), label = round(value, 2)), color = "black")
      }

      p_graphon <- p_graphon +
        ggplot2::scale_y_reverse() +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("Column Blocks") +
        ggplot2::ylab("Row Blocks") +
        ggplot2::ggtitle(x[["net_id"]][net_id]) +
        ggplot2::coord_equal(expand = FALSE)
      return(p_graphon)
    },
    meso = {
      # The below order use mean over all networks to have a consistent display
      if (x$Q[2] == 1) {
        mean_rho <- 1
      } else {
        mean_rho <- matrixStats::rowMeans2(sapply(x$pi, function(pi) pi[[2]]))
      }
      if (is.null(oRow)) {
        oRow <- order(x$alpha %*% mean_rho, decreasing = TRUE)
      }
      #  Once order of tile in block will be fix, need to go back in nullity cond
      if (x$Q[1] == 1) {
        mean_pi <- 1
      } else {
        mean_pi <- matrixStats::rowMeans2(sapply(x$pi, function(pi) pi[[1]]))
      }
      if (is.null(oCol)) {
        oCol <- order(mean_pi %*% x$alpha, decreasing = TRUE)
      }
      p_alpha <- x$alpha[oRow, oCol, drop = FALSE] |>
        t() |>
        reshape2::melt() |>
        ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2("alpha",
          low = "white",
          high = "red",
          limits = c(
            0,
            ifelse(x$distribution == "bernoulli", 1, max(x$alpha))
          )
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "α")) +
        ggplot2::geom_hline(yintercept = seq(x$Q[1]) + .5) +
        ggplot2::geom_vline(xintercept = seq(x$Q[2]) + .5) +
        ggplot2::scale_x_continuous(breaks = seq(x$Q[2])) +
        ggplot2::scale_y_reverse(breaks = seq(x$Q[1])) +
        ggplot2::theme_bw(base_size = 15, base_rect_size = 1, base_line_size = 1) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_fixed(expand = FALSE)

      if (values) {
        p_alpha <- p_alpha +
          ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), color = "black")
      }
      #  scale_y_reverse()
      xl <- ""
      # names(df_pi) <- x$net_id
      if (mixture) {
        df_pi <- purrr::map_dfc(
          seq_along(x$net_id),
          function(m) {
            setNames(
              data.frame(x$pim[[m]][[1]][oRow]),
              x$net_id[m]
            )
          }
        ) |>
          dplyr::mutate(q = seq(x$Q[1])) |>
          tidyr::pivot_longer(cols = -c(q)) |>
          dplyr::mutate(Proportion = value)
        p_pi <-
          df_pi |>
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
            ncol = x$Q[1] %/% 3 + 1,
            byrow = TRUE
          )) +
          ggplot2::ylab("") +
          ggplot2::ylab(xl) +
          ggplot2::xlab("Row proportions") +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90, vjust = .5,
            hjust = 1
          ), aspect.ratio = 1 / x$Q[1])

        df_rho <- purrr::map_dfc(
          seq_along(x$net_id),
          function(m) {
            setNames(
              data.frame(x$pim[[m]][[2]][oCol]),
              x$net_id[m]
            )
          }
        ) |>
          dplyr::mutate(q = seq(x$Q[2])) |>
          tidyr::pivot_longer(cols = -c(q)) |>
          dplyr::mutate(Proportion = value)
        p_rho <- df_rho |>
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
          # Reversing net id to match row
          ggplot2::scale_y_discrete(limits = rev) +
          # Reversing prop order to match the alpha
          ggplot2::scale_x_reverse() +
          ggplot2::guides(fill = ggplot2::guide_legend(
            ncol = x$Q[2] %/% 3 + 1,
            byrow = TRUE
          )) +
          ggplot2::ylab("") +
          ggplot2::ylab(xl) +
          ggplot2::xlab("Column proportions") +
          ggplot2::theme_classic() +
          ggplot2::theme(aspect.ratio = 1 / x$Q[2])
        if (values) {
          p_pi <- p_pi +
            ggplot2::geom_text(ggplot2::aes(label = round(Proportion, 2)),
              position = ggplot2::position_stack(vjust = 0.5),
              color = "black",
              data = subset(df_pi, round(Proportion, 2) > values_min)
            )
          p_rho <- p_rho +
            ggplot2::geom_text(ggplot2::aes(label = round(Proportion, 2)),
              position = ggplot2::position_stack(vjust = 0.5),
              color = "black",
              data = subset(df_rho, round(Proportion, 2) > values_min)
            )
        }
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
      # The below order use net_id parameters
      if (x$Q[2] == 1) {
        mean_rho <- 1
      } else {
        mean_rho <- x[["pi"]][[net_id]][[2]]
      }
      if (is.null(oRow)) {
        oRow <- order(x$alpham[[net_id]] %*% mean_rho, decreasing = TRUE)
      }
      #  Once order of tile in block will be fix, need to go back in nullity cond
      if (x$Q[1] == 1) {
        mean_pi <- 1
      } else {
        mean_pi <- x[["pi"]][[net_id]][[1]]
      }
      if (is.null(oCol)) {
        oCol <- order(mean_pi %*% x$alpham[[net_id]], decreasing = TRUE)
      }

      if (x$Q[1] == 1) {
        ymin <- rep(0, each = x$Q[2])
        ymax <- rep(1, each = x$Q[2])
      } else {
        ymin <- rep(c(0, cumsum(tabulate(x$Z[[net_id]][[1]])[oRow][1:(x$Q[1] - 1)])), each = x$Q[2]) + 0.5
        ymax <- rep(c(cumsum(tabulate(x$Z[[net_id]][[1]])[oRow])), each = x$Q[2]) + 0.5
      }
      if (x$Q[2] == 1) {
        xmin <- rep(0, x$Q[1])
        xmax <- rep(1, x$Q[1])
      } else {
        xmin <- rep(c(0, cumsum(tabulate(x$Z[[net_id]][[2]])[oCol][1:(x$Q[2] - 1)])), x$Q[1]) + 0.5
        xmax <- rep(cumsum(tabulate(x$Z[[net_id]][[2]])[oCol]), x$Q[1]) + 0.5
      }
      if (values) {
        connection_df <- x$alpha[oRow, oCol] |>
          t() |>
          reshape2::melt()
        connection_df <- connection_df |>
          dplyr::arrange(desc(Var2), Var1) |>
          dplyr::mutate(
            xmin = xmin,
            ymin = ymin,
            xmax = xmax,
            ymax = ymax
          )
      }
      Z1_ordered <- ordered(x$Z[[net_id]][[1]], levels = oRow)
      Z2_ordered <- ordered(x$Z[[net_id]][[2]], levels = oCol)
      row_order <- order(Z1_ordered)
      col_order <- order(Z2_ordered)

      block_df <- as.matrix(x$A[[net_id]])[
        row_order,
        col_order
      ] |>
        reshape2::melt()

      p_block <- ggplot2::ggplot(
        data = block_df, ggplot2::aes(
          x = .data$Var2,
          y = rev(.data$Var1),
          fill = .data$value
        )
      ) +
        ggplot2::geom_tile(show.legend = FALSE) +
        # Order will need to reworked to allow to change tile order
        ggplot2::geom_hline(
          yintercept = cumsum(stats::na.omit(tabulate(x$Z[[net_id]][[1]])[rev(oRow)][x$Q[1]:2])) + .5,
          col = "red", linewidth = .5
        ) +
        ggplot2::geom_vline(
          xintercept = cumsum(stats::na.omit(tabulate(x$Z[[net_id]][[2]])[oCol][1:(x$Q[2] - 1)])) + .5,
          col = "red", linewidth = .5
        ) +
        ggplot2::scale_fill_gradient2(high = "black", mid = "white", low = "transparent")
      if (values) {
        p_block <- p_block +
          ggplot2::geom_rect(ggplot2::aes(
            xmin = xmin, ymin = ymin,
            xmax = xmax, ymax = ymax, alpha = value
          ), fill = "red", data = connection_df) +
          ggplot2::guides(alpha = ggplot2::guide_legend(title = "α")) +
          scale_alpha_continuous(limits = c(0, max(x$alpha)))
      }
      p_block <- p_block +
        ggplot2::ylab("") +
        ggplot2::xlab(x$net_id[net_id]) +
        ggplot2::scale_x_discrete(
          breaks = ""
        ) +
        # ggplot2::scale_y_reverse() +
        ggplot2::scale_y_discrete(
          breaks = "",
          limits = rev
        ) +
        ggplot2::coord_equal(expand = FALSE) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(axis.ticks = ggplot2::element_blank())
      return(p_block)
    }
  )
  return(p)
}

#' Plot the state-space exploration plot for a bipartite collection object
#'
#' @param x a bisbmpop object.
#' @param criterion The criterion to plot. Could be "BICL", "ICL" or "vbound".
#' Default is "BICL".
#' @param ... other arguments to pass to the plot.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang arg_match0
#'
#' @return A plot, a ggplot2 object.
#' @export
plot.bisbmpop <- function(x, criterion = "BICL", ...) {
  stopifnot(inherits(x, "bisbmpop"))
  rlang::arg_match0(criterion, c("BICL", "ICL", "vbound"))
  # One value of BIC-L per Q1 Q2
  criterion_df <- x[[criterion]] |>
    reshape2::melt(value.name = "criterion") |>
    dplyr::rename(Q1 = Var1, Q2 = Var2) |>
    dplyr::mutate(Q1 = as.factor(Q1), Q2 = as.factor(Q2)) |>
    # mutate(Q1 = as.factor(Q1), Q2 = as.factor(Q2)) |>
    # Remove -Inf values
    dplyr::filter(criterion > -Inf)

  Qmax <- dim(x$model_list)
  Q1max <- Qmax[1]
  Q2max <- Qmax[2]

  completeness_mat <- outer(seq(1, Q1max), seq(1, Q2max), Vectorize(function(q1, q2) {
    ifelse(is.null(x$model_list[[q1, q2]]$clustering_is_complete), NA, x$model_list[[q1, q2]]$clustering_is_complete)
  }))
  completeness_df <- completeness_mat |>
    reshape2::melt(value.name = "clustering_is_complete") |>
    dplyr::rename(Q1 = Var1, Q2 = Var2) |>
    dplyr::mutate(Q1 = as.factor(Q1), Q2 = as.factor(Q2)) |>
    dplyr::select(Q1, Q2, clustering_is_complete) |>
    dplyr::filter(!is.na(clustering_is_complete)) |>
    dplyr::mutate(clustering_is_complete = ifelse(clustering_is_complete, "Yes", "No"))

  criterion_complete_df <- dplyr::left_join(x = criterion_df, completeness_df, by = c("Q1", "Q2"))
  criterion_complete_df$is_max <- criterion_complete_df$criterion == max(criterion_complete_df$criterion)

  # Convert TRUE FALSE to factor

  p <- ggplot2::ggplot(data = criterion_complete_df) +
    ggplot2::geom_tile(ggplot2::aes(x = Q2, y = Q1, fill = criterion)) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = Q2, y = Q1,
        fill = criterion, color = clustering_is_complete,
        shape = is_max
      ),
      size = 8
    ) +
    ggplot2::scale_shape_manual(name = "Is Max?", values = c(-9, 17)) +
    ggplot2::scale_fill_distiller(name = criterion, palette = "Spectral", direction = 1L) +
    ggplot2::scale_color_brewer(name = "Is Complete?", type = "div") +
    ggplot2::theme(aspect.ratio = Q1max / Q2max) +
    ggplot2::theme_classic()
  p
}

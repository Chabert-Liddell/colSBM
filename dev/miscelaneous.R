binary_list <- function(l) {
#  browser()
  if(length(l) == 1)  return(l)
  # if (length(l) == 2) return(list(l[1],
  #                                 binary_list(l[2])))
  n <- 1 + floor(runif(1, 0, length(l) - 1))
  s <- sample(l, n)
  d <- setdiff(l,s)
#    l[runif(length(l)) > .5]

  if(runif(1) < .3) {
    return (list(l))
  } else {
    return(list(l,
                #list(
                binary_list(s),
                binary_list(d)))#)
  }
}



best_mod <- function(l) {
  if(length(l) == 1)  return(l[[1]])
  if(length(l) == 2) return(l[[2]])
  return(list(best_mod(l[[2]]),
              best_mod(l[[3]])))
}


res <- binary_list(seq(14))
b <- best_mod(res)
b
unlist(b)


res <- estimate_colSBM(foodwebs[4:8], colsbm_model = "pi", directed = TRUE, model = "bernoulli",
                       global_opts = list("verbosity" = 2,"plot_details" = 1, "Q_max" = 10L), nb_run = 1)
res_cl2 <- clusterize_networks(foodwebs, colsbm_model = "iid", directed = TRUE, model = "bernoulli",
                               global_opts = list("verbosity" = 2, Q_max = 10L), nb_run = 1)
lobstr::obj_sizes(res, res_cl)


tibble::tibble(
  Q = seq(length(res$BICL)),
          ICL = res$ICL,
               BICL = res$BICL,
               vbound = res$vbound,
               SBM = res$ICL_sbm) %>%
  tidyr::pivot_longer(cols = -Q, names_to = "Criterion") %>%
  ggplot2::ggplot(ggplot2::aes(x = Q, y = value,
                      linetype = Criterion, color = Criterion,
                      shape = Criterion)) +
  ggplot2::annotate(geom = "rect",
                    xmin = max(which.max(res$BICL) - res$global_opts$depth, 1),
                    xmax = min(which.max(res$BICL) + res$global_opts$depth, length(res$BICL)),
                    ymin = -Inf,
                    ymax = Inf,
                    fill = "gray90"
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 3) +
  ggplot2::ylab("") +
  ggplot2::theme_bw()

library("ggplot2")

filename <- "divergence_modular_to_nested_29-03-23_10:01:08.Rds"
filename <- paste0(getwd(), "/simulation/data/",filename)

current_data <- readRDS(filename)

ggplot(current_data$results) +
    geom_point(aes(x = divergence, y = BICL - sep_LBM_BICL, group = factor(Current_M), col = factor(Current_M))) +
    geom_hline(yintercept = 0)
library("ggplot2")
library("dplyr")


filename <- "divergence_modular_to_nested_with_3_repetitions_29-03-23_13:55:42.Rds"
filename <- paste0(getwd(), "/simulation/data/", filename)

current_data <- readRDS(filename)

current_data$results <- current_data$results %>% mutate(diff_BICL = BICL - sep_LBM_BICL)

current_data$results <- current_data$results %>%
    group_by(Current_M, Network_id, divergence) %>%
    summarize(mean_diff_BICL = mean(diff_BICL), sd_diff_BICL = sd(diff_BICL))

ggplot(current_data$results) +
    aes(x = divergence, y = mean_diff_BICL, group = factor(Current_M), col = factor(Current_M)) +
    geom_point() +
    geom_line() +
    geom_ribbon(
        aes(
            ymin = mean_diff_BICL - sd_diff_BICL,
            ymax = mean_diff_BICL + sd_diff_BICL, 
            fill = factor(Current_M),
            col = NULL
        ),
        alpha = .25
    ) +
    geom_hline(yintercept = 0)
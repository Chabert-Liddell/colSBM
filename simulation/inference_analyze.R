require("ggplot2")
filenames <- list.files(
    path = "./simulation/data/",
    pattern = "inference_testing_2023-05*",
    full.names = TRUE
)
col_id_BICLS <- c(11, 16, 23, 30, 37)
data_list <- lapply(filenames, readRDS)
result_data_frame <- dplyr::bind_rows(data_list)
result_data_frame <- cbind(result_data_frame, preferred_model = sapply(seq_len(nrow(result_data_frame)), function(n) names(which.max(result_data_frame[n, col_id_BICLS]))))

ggplot(data = result_data_frame) +
    aes(x = epsilon_alpha, group = preferred_model, fill = preferred_model) +
    geom_bar()

ggplot(data = result_data_frame) +
    aes(x = , group = preferred_model, fill = preferred_model) +
    geom_bar()


filenames <- list.files(
    path = "./simulation/data/",
    pattern = "model_selection_check*",
    full.names = TRUE
)

data_list <- lapply(filenames, readRDS)
result_data_frame <- dplyr::bind_rows(data_list)

ggplot(data = result_data_frame) +
    aes(x = epsilon_alpha, group = preferred_model, fill = preferred_model) +
    geom_bar()

ggplot(data = result_data_frame) +
    aes(x = , group = preferred_model, fill = preferred_model) +
    geom_bar()


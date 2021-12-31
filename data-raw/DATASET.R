## code to prepare `DATASET` dataset goes here
my_files <- dir(path = "~/Documents/ecological_network/data/food_web/")
my_files <- my_files[1:(length(my_files)-2)]
food_webs <- lapply(
  X = seq_along(my_files),
  FUN = function(i) {
    df <- read.csv(file = paste0("~/Documents/ecological_network/data/food_web/",
                                 my_files[i]), header = TRUE, row.names = 1)
    A <- as.matrix(df)
    return(list(
      net = A,
      nr = nrow(A),
      nc = ncol(A),
      dens = mean(A),
      id = stringr::str_sub(my_files[i], 1, -5))
    )
  }
)




sub_fw <- purrr::map(food_webs, "net")[11:18]
sub_id <- purrr::map_chr(food_webs, "id")[11:18]
site_names <- c("M_Martins", "NC_Cooper", "NC_Herlzier", "NZ_Venlaw",
                "NZ_Berwick", "NZ_North_Col", "NZ_Powder", "NZ_Trib_C" )
foodwebs <- sub_fw
names(foodwebs) <- site_names
usethis::use_data(foodwebs, overwrite = TRUE)

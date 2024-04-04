future_is_installed <- ("future.apply" %in% installed.packages())

os_is_windows <- (toupper(Sys.info()["sysname"]) == "WINDOWS")

bettermc_is_installed <- ("bettermc" %in% installed.packages())

parallel_is_installed <- ("parallel" %in% installed.packages())

M <- 5
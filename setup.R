pkgs <- c("dplyr", "tidyr", "purrr", "readr", "readxl", "janitor", "httr2", "haven",
          "countrycode", "WDI", "fixest")

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv/activate.R")) renv::init(bare = TRUE)
renv::install(pkgs)
renv::snapshot(packages = pkgs, prompt = FALSE)

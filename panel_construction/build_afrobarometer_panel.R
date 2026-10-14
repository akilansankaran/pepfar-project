suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(countrycode)
})

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
final_dir <- "data/final"

afro <- read_csv(file.path(raw_dir, "afrobarometer_harmonised.csv"), show_col_types = FALSE) |>
  mutate(iso3 = countrycode(country, "country.name", "iso3c", warn = FALSE))

country_panel <- read_csv(file.path(final_dir, "country_panel.csv"), show_col_types = FALSE) |>
  select(iso3, year, PEPFAR, PEPFAR_focus, pepfar_first_year, hiv_prev_baseline,
         pepfar_per_capita, gdppc_growth)

# Respondent-level file keeps the survey weights for individual models; the country-round
# file collapses to weighted means for comparison with the macro panel.
respondents <- afro |>
  left_join(country_panel, by = c("iso3", "survey_year" = "year")) |>
  mutate(post = as.integer(survey_year >= 2004))

country_round <- respondents |>
  group_by(iso3, round, survey_year, PEPFAR, PEPFAR_focus, post) |>
  summarise(across(starts_with("econ_"), ~ weighted.mean(.x, weight, na.rm = TRUE)),
            n_resp = n(), .groups = "drop")

dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(respondents, file.path(final_dir, "afrobarometer_respondents.csv"))
write_csv(country_round, file.path(final_dir, "afrobarometer_country_round.csv"))

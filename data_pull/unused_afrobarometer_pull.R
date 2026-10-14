suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(haven)
})

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
afro_dir <- file.path(raw_dir, "afrobarometer")

# Merged-round SPSS files from afrobarometer.org (free, registration required), saved as
# r1.sav ... r9.sav. Round 1 (1999-2001) and round 2 (2002-2003) are the only pre-PEPFAR
# waves, so they anchor any before/after comparison.
round_years <- c(r1 = 2000, r2 = 2003, r3 = 2005, r4 = 2008, r5 = 2012,
                 r6 = 2015, r7 = 2018, r8 = 2021, r9 = 2023)

# Question numbering changes across rounds; each entry maps a harmonised name to the
# variable holding it in that round, and must be checked against each round's codebook
# before use. Rounds missing an item are left NA.
item_map <- list(
  econ_country_now   = c(r1 = "q6",  r2 = "q4a", r3 = "q4a", r4 = "q4a", r5 = "q3a",
                         r6 = "q4a", r7 = "q4a", r8 = "q4a", r9 = "q4a"),
  econ_country_ahead = c(r2 = "q7a", r3 = "q6a", r4 = "q6a", r5 = "q6a",
                         r6 = "q6a", r7 = "q6a", r8 = "q6a", r9 = "q6a"),
  econ_own_now       = c(r1 = "q7",  r2 = "q4b", r3 = "q4b", r4 = "q4b", r5 = "q3b",
                         r6 = "q4b", r7 = "q4b", r8 = "q4b", r9 = "q4b"),
  econ_own_ahead     = c(r2 = "q7b", r3 = "q6b", r4 = "q6b", r5 = "q6b",
                         r6 = "q6b", r7 = "q6b", r8 = "q6b", r9 = "q6b")
)

read_round <- function(round) {
  path <- file.path(afro_dir, paste0(round, ".sav"))
  if (!file.exists(path)) {
    warning("missing ", path)
    return(NULL)
  }
  d <- read_sav(path)
  names(d) <- tolower(names(d))
  out <- tibble(
    round = round,
    survey_year = round_years[[round]],
    country = as_factor(d$country),
    region = if ("region" %in% names(d)) as.character(as_factor(d$region)) else NA_character_,
    weight = if ("withinwt" %in% names(d)) as.numeric(d$withinwt) else 1
  )
  for (item in names(item_map)) {
    var <- item_map[[item]][round]
    # codes above 5 are don't know, refused, and not applicable
    out[[item]] <- if (!is.na(var) && var %in% names(d)) {
      v <- as.numeric(zap_labels(d[[var]]))
      replace(v, v > 5 | v < 1, NA)
    } else NA_real_
  }
  out
}

afro <- map(names(round_years), read_round) |> compact() |> bind_rows()
write_csv(afro, file.path(raw_dir, "afrobarometer_harmonised.csv"))

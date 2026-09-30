suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(countrycode)
})

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
final_dir <- "data/final"

program_start <- 2004
baseline_year <- 2003

# The fifteen original focus countries named under the 2003 Leadership Act.
focus_countries <- c("BWA", "CIV", "ETH", "GUY", "HTI", "KEN", "MOZ", "NAM",
                     "NGA", "RWA", "ZAF", "TZA", "UGA", "VNM", "ZMB")

read_raw <- function(name) {
  path <- file.path(raw_dir, name)
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}

wdi    <- read_raw("wdi.csv")
fa     <- read_raw("foreignassistance_hiv.csv")
crs    <- read_raw("oecd_crs_hiv.csv")
wpp    <- read_raw("un_wpp.csv")
income <- read_raw("wb_income_history.csv")
gbd    <- read_raw("ihme_gbd.csv")

if (!is.null(gbd)) {
  gbd <- gbd |>
    mutate(iso3 = countrycode(location_name, "country.name", "iso3c", warn = FALSE)) |>
    filter(!is.na(iso3)) |>
    select(-location_name)
}

panel <- wdi |>
  left_join(fa, by = c("iso3", "year")) |>
  left_join(crs, by = c("iso3", "year")) |>
  left_join(wpp, by = c("iso3", "year")) |>
  left_join(income, by = c("iso3", "year"))
if (!is.null(gbd)) panel <- left_join(panel, gbd, by = c("iso3", "year"))

# The sample is every country classified low or middle income at baseline with a usable
# growth series; this is the rule that should reproduce the 157-country panel, and the
# resulting list is checked against the replication materials in analysis/replication.
baseline_income <- income |>
  filter(year == baseline_year) |>
  select(iso3, income_baseline = income_group)

sample_iso3 <- panel |>
  left_join(baseline_income, by = "iso3") |>
  filter(income_baseline %in% c("L", "LM", "UM")) |>
  group_by(iso3) |>
  filter(sum(!is.na(gdppc_growth)) >= 10) |>
  distinct(iso3) |>
  pull(iso3)

panel <- panel |>
  filter(iso3 %in% sample_iso3) |>
  left_join(baseline_income, by = "iso3") |>
  mutate(pepfar_disb_const = coalesce(pepfar_disb_const, 0))

treatment <- panel |>
  group_by(iso3) |>
  summarise(
    cum_pepfar = sum(pepfar_disb_const[year >= program_start]),
    pepfar_first_year = suppressWarnings(min(year[pepfar_disb_const > 0 & year >= program_start])),
    .groups = "drop"
  ) |>
  mutate(pepfar_first_year = ifelse(is.finite(pepfar_first_year), pepfar_first_year, NA))

baseline_burden <- panel |>
  filter(year == baseline_year) |>
  select(iso3, hiv_prev_baseline = hiv_prev_15_49,
         any_of(c(hiv_daly_baseline = "hiv_dalys", hiv_death_baseline = "hiv_deaths")))

panel <- panel |>
  left_join(treatment, by = "iso3") |>
  left_join(baseline_burden, by = "iso3") |>
  mutate(
    PEPFAR = as.integer(cum_pepfar > 0),
    PEPFAR_focus = as.integer(iso3 %in% focus_countries),
    Time = as.integer(year >= program_start),
    Interaction = PEPFAR * Time,
    Interaction_focus = PEPFAR_focus * Time,
    event_time = year - pepfar_first_year,
    pepfar_per_capita = pepfar_disb_const / pop_total,
    region_wb = region,
    country = countrycode(iso3, "iso3c", "country.name")
  ) |>
  select(-region) |>
  arrange(iso3, year)

message(n_distinct(panel$iso3), " countries, ",
        sum(panel$PEPFAR[!duplicated(panel$iso3)]), " PEPFAR recipients")

dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(panel, file.path(final_dir, "country_panel.csv"))

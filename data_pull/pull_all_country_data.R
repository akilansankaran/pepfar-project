suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(janitor)
  library(WDI)
  library(httr2)
})

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

start_year <- 1990
end_year <- 2023

download_once <- function(url, dest) {
  if (!file.exists(dest)) {
    message("downloading ", basename(dest))
    request(url) |>
      req_retry(max_tries = 3) |>
      req_perform(path = dest)
  }
  dest
}

# Outcomes, the controls used in the replication, and the capital-inflow measures used to
# address reverse causality in the disease-burden extension.
wdi_indicators <- c(
  gdp_growth        = "NY.GDP.MKTP.KD.ZG",
  gdppc_growth      = "NY.GDP.PCAP.KD.ZG",
  gdppc_const       = "NY.GDP.PCAP.KD",
  life_exp          = "SP.DYN.LE00.IN",
  crude_death_rate  = "SP.DYN.CDRT.IN",
  hiv_prev_15_49    = "SH.DYN.AIDS.ZS",
  gross_cap_form    = "NE.GDI.TOTL.ZS",
  trade_open        = "NE.TRD.GNFS.ZS",
  inflation         = "FP.CPI.TOTL.ZG",
  pop_growth        = "SP.POP.GROW",
  pop_total         = "SP.POP.TOTL",
  primary_enroll    = "SE.PRM.ENRR",
  fdi_inflow_gdp    = "BX.KLT.DINV.WD.GD.ZS",
  remittances_gdp   = "BX.TRF.PWKR.DT.GD.ZS",
  oda_gni           = "DT.ODA.ODAT.GN.ZS",
  ext_debt_gni      = "DT.DOD.DECT.GN.ZS"
)

pull_wdi <- function() {
  WDI(country = "all", indicator = wdi_indicators,
      start = start_year, end = end_year, extra = TRUE) |>
    as_tibble() |>
    filter(region != "Aggregates") |>
    select(iso3 = iso3c, year, region, income_wdi = income, lending,
           all_of(names(wdi_indicators))) |>
    write_csv(file.path(raw_dir, "wdi.csv"))
}

# The complete foreignassistance.gov extract is large (several GB); PEPFAR spending is
# identified through the HIV/AIDS sector, which is where State/PEPFAR obligations are coded.
fa_url <- Sys.getenv(
  "FA_GOV_URL",
  "https://s3.amazonaws.com/files.explorer.devtechlab.com/us_foreign_aid_complete.csv"
)

pull_foreignassistance <- function() {
  path <- download_once(fa_url, file.path(raw_dir, "us_foreign_aid_complete.csv"))
  read_csv(path, show_col_types = FALSE, guess_max = 1e5) |>
    clean_names() |>
    filter(us_sector_name == "HIV/AIDS",
           transaction_type_name == "Disbursements") |>
    group_by(iso3 = country_code, year = fiscal_year) |>
    summarise(pepfar_disb_const = sum(constant_dollar_amount, na.rm = TRUE),
              .groups = "drop") |>
    mutate(year = as.integer(year)) |>
    write_csv(file.path(raw_dir, "foreignassistance_hiv.csv"))
}

# CRS through the OECD SDMX endpoint. Purpose code 13040 is STD control including HIV/AIDS;
# all donors are kept so that non-US HIV aid can enter as a competing intervention.
# The key must follow the dimension order of DSD_CRS; check it against the dataflow
# definition if OECD revises the structure.
crs_url <- paste0(
  "https://sdmx.oecd.org/public/rest/data/OECD.DCD.FSD,DSD_CRS@DF_CRS,1.4/",
  "..13040.100._T._T.D.Q._T..",
  "?startPeriod=", start_year, "&endPeriod=", end_year,
  "&dimensionAtObservation=AllDimensions&format=csvfilewithlabels"
)

pull_oecd_crs <- function() {
  path <- download_once(crs_url, file.path(raw_dir, "oecd_crs_13040.csv"))
  read_csv(path, show_col_types = FALSE) |>
    clean_names() |>
    transmute(donor, iso3 = recipient, year = as.integer(time_period),
              hiv_aid_disb = obs_value) |>
    group_by(iso3, year) |>
    summarise(hiv_aid_all_donors = sum(hiv_aid_disb, na.rm = TRUE),
              hiv_aid_us = sum(hiv_aid_disb[donor == "USA"], na.rm = TRUE),
              .groups = "drop") |>
    mutate(hiv_aid_non_us = hiv_aid_all_donors - hiv_aid_us) |>
    write_csv(file.path(raw_dir, "oecd_crs_hiv.csv"))
}

# GBD has no open bulk API; exports from the GBD Results Tool are placed in raw_dir/ihme_gbd/
# by hand. Expected query: cause HIV/AIDS and all causes, measures deaths, DALYs, prevalence
# and incidence, metric rate, age-standardised, both sexes, all countries and years.
pull_ihme_gbd <- function() {
  files <- list.files(file.path(raw_dir, "ihme_gbd"), pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    warning("no GBD exports in ", file.path(raw_dir, "ihme_gbd"), "; skipping")
    return(invisible(NULL))
  }
  lapply(files, read_csv, show_col_types = FALSE) |>
    bind_rows() |>
    filter(sex_name == "Both", metric_name == "Rate",
           age_name == "Age-standardized") |>
    mutate(var = paste(
      ifelse(cause_name == "All causes", "allcause", "hiv"),
      tolower(gsub("[^A-Za-z]+.*", "", measure_name)), sep = "_"
    )) |>
    select(location_name, year, var, val) |>
    pivot_wider(names_from = var, values_from = val) |>
    write_csv(file.path(raw_dir, "ihme_gbd.csv"))
}

wpp_url <- paste0(
  "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/",
  "WPP2024_Demographic_Indicators_Medium.csv.gz"
)

pull_un_population <- function() {
  path <- download_once(wpp_url, file.path(raw_dir, "wpp2024_demographic.csv.gz"))
  read_csv(path, show_col_types = FALSE, guess_max = 1e5) |>
    filter(!is.na(ISO3_code), Variant == "Medium",
           between(Time, start_year, end_year)) |>
    transmute(iso3 = ISO3_code, year = as.integer(Time),
              un_pop_thousands = TPopulation1July,
              un_pop_density = PopDensity,
              un_median_age = MedianAgePop,
              un_life_exp = LEx,
              un_under5_mort = Q5,
              un_net_migration = NetMigrations) |>
    write_csv(file.path(raw_dir, "un_wpp.csv"))
}

# Historical (year-by-year) income groups, needed because treatment and control pools are
# drawn by classification at baseline rather than by today's grouping.
oghist_url <- Sys.getenv(
  "WB_OGHIST_URL",
  "https://databank.worldbank.org/data/download/site-content/OGHIST.xlsx"
)

pull_wb_classification <- function() {
  path <- download_once(oghist_url, file.path(raw_dir, "wb_oghist.xlsx"))
  hist <- read_excel(path, sheet = "Country Analytical History", skip = 5,
                     col_names = FALSE)
  years <- as.integer(unlist(read_excel(path, sheet = "Country Analytical History",
                                        range = "C6:AZ6", col_names = FALSE)))
  hist <- hist[, seq_len(2 + sum(!is.na(years)))]
  names(hist) <- c("iso3", "country", paste0("y", years[!is.na(years)]))
  hist |>
    filter(nchar(iso3) == 3) |>
    pivot_longer(starts_with("y"), names_to = "year", values_to = "income_group") |>
    mutate(year = as.integer(sub("y", "", year)),
           income_group = na_if(trimws(income_group), "..")) |>
    select(-country) |>
    write_csv(file.path(raw_dir, "wb_income_history.csv"))
}

pull_wdi()
pull_foreignassistance()
pull_oecd_crs()
pull_ihme_gbd()
pull_un_population()
pull_wb_classification()

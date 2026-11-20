# Step 5 inputs: GNI per capita (PPP) and the US GDP deflator from the World Bank WDI
# Author: Rishad Vaghaiwalla
#
# Used by analysis/parallel_trends/cost_benefit_income_groups.R for two things:
#   1. Converting the World Bank income thresholds (current US$ of the year they are based
#      on) into the panel's constant 2022 US$, with the same US GDP deflator as
#      build_pepfar_panel_FIXED.R: value_2022 = value_t * (deflator_2022 / deflator_t).
#   2. The VSL cross-check (Robinson, Hammitt and O'Keeffe 2019), which extrapolates from the
#      US using GNI per capita at purchasing power parity:
#        gni_pc_ppp   GNI per capita, PPP (current international $)   NY.GNP.PCAP.PP.CD
#
# Same download helper as data_pull/pull_wdi_mechanisms.R (one indicator at a time, with
# retries, regional aggregates dropped).
#
# Run from the repository root:  Rscript data_pull/pull_wdi_valuation.R
# Writes data/intermediate/wdi_valuation_2000_2022.csv and data/intermediate/us_gdp_deflator.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(WDI)
})

YEARS      <- 2000:2022
DEFL_YEARS <- 2000:2025     # the deflator also needs the threshold year (FY2027 uses 2025 GNI)
out_file   <- "data/intermediate/wdi_valuation_2000_2022.csv"
defl_file  <- "data/intermediate/us_gdp_deflator.csv"

wdi_one <- function(name, code, country = "all", years = YEARS, tries = 4) {
  for (i in seq_len(tries)) {
    d <- tryCatch(WDI(country = country, indicator = setNames(code, name),
                      start = min(years), end = max(years)),
                  error = function(e) NULL)
    if (!is.null(d) && name %in% names(d)) {
      return(d[!is.na(d$iso3c) & grepl("^[A-Z]{3}$", d$iso3c), c("iso3c", "year", name)])
    }
    message("  retry ", code, " (", i, ")"); Sys.sleep(5 * i)
  }
  stop("WDI download failed for ", code)
}

gni  <- wdi_one("gni_pc_ppp", "NY.GNP.PCAP.PP.CD")
defl <- wdi_one("deflator", "NY.GDP.DEFL.ZS", country = "US", years = DEFL_YEARS) |>
  select(year, deflator) |> filter(!is.na(deflator)) |> arrange(year)

write_csv(gni, out_file)
write_csv(defl, defl_file)
cat("Wrote", out_file, "(", nrow(gni), "rows ) and", defl_file, "\n")
cat("US deflator available through", max(defl$year), "\n")
cat("Share of country-years with GNI per capita (PPP):", round(mean(!is.na(gni$gni_pc_ppp)), 3), "\n")

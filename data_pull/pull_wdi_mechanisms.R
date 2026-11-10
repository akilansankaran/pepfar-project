# Investment and savings mechanism variables from the World Bank WDI
# Author: Anthony Babu
#
# Mechanism being tested: PEPFAR lowers adult mortality -> longer expected working lives
# -> more saving and investment (domestic and foreign). These series are the outcomes
# for analysis/mechanisms/fdi_panel.R. They are not in pepfar_panel_1990_2022_FIXED.csv,
# so this script pulls them separately and writes a file that merges on iso3c + year.
#
#   fdi_pct_gdp     FDI, net inflows (% of GDP)                    BX.KLT.DINV.WD.GD.ZS
#   fdi_pc          FDI, net inflows per person, constant 2022 US$ BX.KLT.DINV.CD.WD / pop_total
#   savings_pct_gdp Gross domestic savings (% of GDP)              NY.GDS.TOTL.ZS
#   gfcf_pct_gdp    Gross fixed capital formation (% of GDP)       NE.GDI.FTOT.ZS
#   gross_sav_pc    Gross savings per person, constant 2022 US$     NY.GNS.ICTR.CD / pop_total
#   gfcf_pc         Gross fixed capital formation per person,
#                   constant 2022 US$                              NE.GDI.FTOT.CD / pop_total
#
# The two per-person series are the "not normed by GDP" versions: if PEPFAR raised GDP, a
# share of GDP can stay flat while saving or investment rises. They use the current-US$
# series (not local currency, .CN) so every country is in the same dollars and gets the
# same deflator as the rest of the panel.
#
# fdi_pc, gross_sav_pc and gfcf_pc use the same conventions as build_pepfar_panel_FIXED.R: current US$ deflated
# with the US GDP deflator (NY.GDP.DEFL.ZS), value_constant = value_current *
# (deflator_2022 / deflator_t), and divided by the panel's UN WPP population.
#
# Run from the repository root:  Rscript data_pull/pull_wdi_mechanisms.R
# Writes data/intermediate/wdi_mechanisms_1990_2022.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(WDI)
})

YEARS     <- 1990:2022
BASE_YEAR <- 2022
out_file  <- "data/intermediate/wdi_mechanisms_1990_2022.csv"

panel <- read_csv("data/intermediate/pepfar_panel_1990_2022_FIXED.csv",
                  col_select = c(iso3c, year, pop_total), show_col_types = FALSE)

codes <- c(fdi_pct_gdp     = "BX.KLT.DINV.WD.GD.ZS",
           fdi_cur         = "BX.KLT.DINV.CD.WD",
           savings_pct_gdp = "NY.GDS.TOTL.ZS",
           gfcf_pct_gdp    = "NE.GDI.FTOT.ZS",
           gross_sav_cur   = "NY.GNS.ICTR.CD",
           gfcf_cur        = "NE.GDI.FTOT.CD")

# One indicator at a time with retries; the World Bank API is often slow
wdi_one <- function(name, code, country = "all", tries = 4) {
  for (i in seq_len(tries)) {
    d <- tryCatch(WDI(country = country, indicator = setNames(code, name),
                      start = min(YEARS), end = max(YEARS)),
                  error = function(e) NULL)
    if (!is.null(d) && name %in% names(d)) {
      # drop regional aggregates (no ISO3 code), which otherwise duplicate rows
      return(d[!is.na(d$iso3c) & grepl("^[A-Z]{3}$", d$iso3c), c("iso3c", "year", name)])
    }
    message("  retry ", code, " (", i, ")"); Sys.sleep(5 * i)
  }
  stop("WDI download failed for ", code)
}

wdi <- Reduce(\(a, b) full_join(a, b, by = c("iso3c", "year")),
              Map(wdi_one, names(codes), codes))

defl <- wdi_one("deflator", "NY.GDP.DEFL.ZS", country = "US") |> select(year, deflator)
d_base <- defl$deflator[defl$year == BASE_YEAR]

out <- panel |>
  left_join(wdi, by = c("iso3c", "year")) |>
  left_join(defl, by = "year") |>
  mutate(across(c(fdi_cur, gross_sav_cur, gfcf_cur), \(x) x * (d_base / deflator) / pop_total,
                .names = "{sub('_cur', '_pc', .col)}")) |>
  select(iso3c, year, fdi_pct_gdp, fdi_pc, savings_pct_gdp, gfcf_pct_gdp, gross_sav_pc, gfcf_pc)

stopifnot(nrow(out) == nrow(panel))  # one row per panel country-year
write_csv(out, out_file)
cat("Wrote", out_file, "\n")
print(summarise(out, across(-c(iso3c, year), \(x) mean(!is.na(x)))))

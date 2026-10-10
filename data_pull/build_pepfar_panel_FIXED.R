# =============================================================================
# build_pepfar_panel.R
#
# Builds a country-year panel (157 countries x 1990-2022) for the PEPFAR
# economic-growth (Crown) and mortality (Kates) analyses.
#
# Sources
#   - World Bank WDI (API, via the WDI package)
#   - foreignassistance.gov complete dataset (PEPFAR / US HIV funding)
#   - OECD CRS (SDMX API) - donor health ODA, non-PEPFAR
#   - GBD results tool (manual download) - HIV prevalence gap-fill
#   - UN WPP - total population
#   - World Bank OGHIST - historical income classification
#
# Conventions
#   - All dollar values converted to constant 2022 USD using the US GDP
#     deflator (WDI NY.GDP.DEFL.ZS, USA):
#         value_constant = value_current * (deflator_2022 / deflator_t)
#     Every series is pulled in CURRENT dollars and deflated here, so all money
#     variables share one deflator and one base year.
#   - PEPFAR amounts are by US fiscal year and exclude COVID-19 funding.
#   - Maternal mortality only exists from 2000 onward (set NA before 2000).
#
# Output (in OUT_DIR)
#   pepfar_panel_1990_2022.csv / .rds   - the panel
#   codebook.csv                        - variable list with source
#   sample_countries.csv                - the country list actually used
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(WDI)
  library(countrycode)
})

# ---- Configuration ----------------------------------------------------------

YEARS        <- 1990:2022
BASE_YEAR    <- 2022           # constant-dollar base year
BASELINE_YR  <- 2003           # pre-PEPFAR baseline year for baseline covariates
PEPFAR_START <- 2004           # first PEPFAR fiscal year
N_EXPECTED   <- 157

RAW_DIR <- "data/raw"
OUT_DIR <- "data/clean"
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Study country list (157 countries) with treatment groups:
#   COP-PEPFAR (31), Other PEPFAR recipient (59), Control (67)
SAMPLE_FILE <- file.path(RAW_DIR, "pepfar_country_list.xlsx")
# The list carries Serbia (SRB), Montenegro (MNE) and the former Serbia and
# Montenegro (SCG) as separate units. SCG only has data where a source still
# reports under that code (mainly pre-2006 aid flows); WDI/WPP have none.
#
# RESOLVED (checked [DATE]): the paper's control threshold is "<$1M OR
# <$0.05 per capita" cumulative PEPFAR-linked disbursement -- an OR, not a
# strict $1M cutoff. Re-checking the four countries an earlier, incomplete
# $1M-only pass had flagged, against BOTH conditions, for both windows:
#   - VEN: $0 / $0.00pc through 2018; $7.01M / $0.249pc through 2022
#     -> valid Control for Crown (2004-2018), NOT valid for Kates (2004-2022).
#   - COL: $1.34M / $0.027pc through 2018 (per-capita saves it);
#     $4.50M / $0.087pc through 2022 (fails both)
#     -> valid Control for Crown, NOT valid for Kates.
#   - PAN: $120K / $0.029pc through 2018; $4.45M / $1.010pc through 2022
#     (small population, so per-capita spikes hard)
#     -> valid Control for Crown, NOT valid for Kates.
#   - EGY: $3.56M / $0.034pc through 2018; $3.56M / $0.032pc through 2022
#     (per-capita stays under threshold in BOTH windows despite the total
#     exceeding $1M) -> valid Control for Crown AND Kates.
# VEN/COL/PAN's violation is driven by 2019-2022 Dept. of State Global
# Health Programs funding tied to PEPFAR's regional response to the
# Venezuelan migrant/refugee crisis (Colombia and Panama are the two
# largest host countries) -- this money simply didn't exist as of Crown's
# 2004-2018 window.
# ACTION TAKEN: added grp_control_kates (0/1) below, identical to
# grp_control except VEN/COL/PAN are set to 0. Use grp_control_kates (not
# grp_control) when restricting the sample for any Kates-window (2004-2022)
# analysis; continue using grp_control for Crown-window (2004-2018) work.
# NOTE: only these four were checked this way -- the other 63 "Control"
# countries have not been individually re-verified against the Kates-window
# threshold. Provenance of the Group column itself is still undocumented.

# foreignassistance.gov "complete" dataset (several GB). Download once.
FA_URL  <- "https://s3.amazonaws.com/files.explorer.devtechlab.com/us_foreign_aid_complete.csv"
FA_FILE <- file.path(RAW_DIR, "us_foreign_aid_complete.csv")
FA_HIV_FILE <- file.path(RAW_DIR, "fa_hiv_disbursements.csv")   # filtered cache

# UN WPP 2024 demographic indicators (medium variant)
WPP_URL  <- "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/WPP2024_Demographic_Indicators_Medium.csv.gz"
WPP_FILE <- file.path(RAW_DIR, "WPP2024_Demographic_Indicators_Medium.csv.gz")

# World Bank historical income classifications
OGHIST_URL  <- "https://datacatalogfiles.worldbank.org/ddh-published/0037712/DR0090754/OGHIST.xlsx"
OGHIST_FILE <- file.path(RAW_DIR, "OGHIST.xlsx")

# GBD HIV prevalence, downloaded by hand from https://vizhub.healthdata.org/gbd-results/
#   Measure: Prevalence | Metric: Percent | Cause: HIV/AIDS | Age: 15-49 years
#   Sex: Both | Location: all countries | Year: 1990-2022
# Save as CSV here. Skipped (with a warning) if not present.
GBD_FILE <- file.path(RAW_DIR, "gbd_hiv_prevalence_15_49.csv")

# OECD CRS SDMX endpoint. Key order for DSD_CRS@DF_CRS:
#   DONOR.RECIPIENT.SECTOR.MEASURE.CHANNEL.MODALITY.FLOW_TYPE.PRICE_BASE.MD_DIM.MD_ID.UNIT_MEASURE
# PRICE_BASE V = current prices (deflated below); FLOW_TYPE D = disbursements.
# CRS_DONOR: "DAC" = DAC member total (as in the reference query). To include
# multilaterals, swap in the all-donors aggregate code from the Data Explorer.
CRS_BASE  <- "https://sdmx.oecd.org/dcd-public/rest/data/OECD.DCD.FSD,DSD_CRS@DF_CRS,"
CRS_DONOR <- "DAC"

# ---- Helpers ----------------------------------------------------------------

download_once <- function(url, dest, mode = "wb") {
  if (!file.exists(dest)) {
    message("Downloading ", basename(dest), " ...")
    options(timeout = max(3600, getOption("timeout")))
    download.file(url, dest, mode = mode, quiet = FALSE)
  }
  invisible(dest)
}

clean_names <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", x)
}

to_iso3 <- function(name) {
  suppressWarnings(countrycode(name, "country.name", "iso3c",
                               custom_match = c("Kosovo" = "XKX",
                                                "Micronesia" = "FSM")))
}

deflate <- function(value_current, year, defl) {
  # value_constant = value_current * (deflator_base_year / deflator_current_year)
  d_t    <- defl$deflator[match(year, defl$year)]
  d_base <- defl$deflator[defl$year == BASE_YEAR]
  value_current * (d_base / d_t)
}

# WDI pull, one indicator at a time with retries; cached to RAW_DIR so reruns
# don't hit the (often slow) World Bank API.
wdi_pull <- function(codes, country = "all", tries = 4) {
  cache <- file.path(RAW_DIR, paste0("wdi_", country, ".rds"))
  have  <- if (file.exists(cache)) readRDS(cache) else NULL
  todo  <- setdiff(names(codes), names(have))
  for (v in todo) {
    for (i in seq_len(tries)) {
      d <- tryCatch(WDI(country = country, indicator = codes[v],
                        start = min(YEARS), end = max(YEARS), extra = TRUE),
                    error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(d) && v %in% names(d)) break
      message("  retry ", v, " (", i, ")"); Sys.sleep(5 * i)
    }
    if (is.null(d) || !v %in% names(d)) stop("WDI download failed for ", codes[[v]])
    d <- d[!is.na(d$iso3c) & grepl("^[A-Z]{3}$", d$iso3c) &
             !d$region %in% "Aggregates", c("iso3c", "year", v)]
    have <- if (is.null(have)) d else full_join(have, d, by = c("iso3c", "year"))
    saveRDS(have, cache)
  }
  have
}

safe_log1p <- function(x) {
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(x) & x > -1
  out[ok] <- log1p(x[ok])
  out
}

# =============================================================================
# 1. US GDP deflator (for constant 2022 USD)
# =============================================================================

defl <- wdi_pull(c(deflator = "NY.GDP.DEFL.ZS"), country = "US") |>
  select(year, deflator) |>
  arrange(year)
stopifnot(!anyNA(defl$deflator), BASE_YEAR %in% defl$year)

# =============================================================================
# 2. World Development Indicators
# =============================================================================

wdi_codes <- c(
  # economic
  gdp_pc_cur        = "NY.GDP.PCAP.CD",      # GDP per capita, current US$
  gdp_cur           = "NY.GDP.MKTP.CD",      # GDP, current US$
  gdp_pc_growth     = "NY.GDP.PCAP.KD.ZG",   # GDP per capita growth, annual %
  gdp_growth        = "NY.GDP.MKTP.KD.ZG",   # GDP growth, annual %
  emp_female        = "SL.EMP.TOTL.SP.FE.ZS",# employment-to-pop, 15+, female (modeled ILO)
  emp_male          = "SL.EMP.TOTL.SP.MA.ZS",# employment-to-pop, 15+, male (modeled ILO)
  oos_prim_female   = "SE.PRM.UNER.FE",      # children out of school, primary, female (count)
  oos_prim_male     = "SE.PRM.UNER.MA",      # children out of school, primary, male (count)
  pop_prim_female   = "SP.PRM.TOTL.FE.IN",   # primary school-age population, female
  pop_prim_male     = "SP.PRM.TOTL.MA.IN",   # primary school-age population, male
  # mortality
  mort_allcause     = "SP.DYN.CDRT.IN",      # crude death rate per 1,000 people
  mort_child        = "SH.DYN.MORT",         # under-5 mortality per 1,000 live births
  mort_maternal     = "SH.STA.MMRT",         # maternal mortality ratio per 100,000 live births
  # baseline covariates
  life_exp          = "SP.DYN.LE00.IN",
  fertility         = "SP.DYN.TFRT.IN",
  urban_pct         = "SP.URB.TOTL.IN.ZS",
  sec_enroll_gross  = "SE.SEC.ENRR",
  hiv_prev_wdi      = "SH.DYN.AIDS.ZS",      # % of population ages 15-49
  dom_hexp_gov_ppp  = "SH.XPD.GHED.PP.CD",   # domestic general govt health exp per capita, PPP current intl $
  dom_hexp_pvt_ppp  = "SH.XPD.PVTD.PP.CD"    # domestic private health exp per capita, PPP current intl $
)

wdi_raw <- wdi_pull(wdi_codes)

wdi <- wdi_raw |>
  select(iso3c, year, all_of(names(wdi_codes))) |>
  as_tibble()

# =============================================================================
# 3. UN WPP total population
# =============================================================================

download_once(WPP_URL, WPP_FILE)
wpp <- read_csv(WPP_FILE, col_select = c(ISO3_code, Time, Variant, TPopulation1July),
                show_col_types = FALSE) |>
  filter(!is.na(ISO3_code), Time %in% YEARS, Variant == "Medium") |>
  transmute(iso3c = ISO3_code, year = as.integer(Time),
            pop_total = TPopulation1July * 1000) |>   # WPP reports thousands
  distinct(iso3c, year, .keep_all = TRUE)

# =============================================================================
# 4. World Bank historical income classification (OGHIST)
# =============================================================================

download_once(OGHIST_URL, OGHIST_FILE)
og <- as.data.frame(suppressMessages(
  read_excel(OGHIST_FILE, sheet = "Country Analytical History", col_names = FALSE)))

# The "Data for calendar year :" row holds the data years; country rows have
# an ISO3 code in column 1. Classes: L, LM, UM, H (".." = not classified).
yr_row  <- which(apply(og, 1, \(r) any(grepl("calendar year", r, ignore.case = TRUE))))[1]
yr_vals <- suppressWarnings(as.integer(unlist(og[yr_row, ])))
yr_cols <- which(!is.na(yr_vals))

og_c <- og[grepl("^[A-Z]{3}$", og[[1]]), c(1, yr_cols)]
names(og_c) <- c("iso3c", yr_vals[yr_cols])

income_class <- og_c |>
  pivot_longer(-iso3c, names_to = "year", values_to = "wb_income") |>
  mutate(year = as.integer(year),
         wb_income = toupper(trimws(gsub("\\*", "", wb_income))),
         wb_income = if_else(wb_income %in% c("L", "LM", "UM", "H"), wb_income, NA_character_)) |>
  filter(year %in% YEARS) |>
  mutate(wb_income = factor(wb_income, levels = c("L", "LM", "UM", "H")))

# =============================================================================
# 5. foreignassistance.gov: PEPFAR spending + pre-2004 US HIV funding dummy
# =============================================================================

# The complete file is ~3.7 GB. DuckDB scans it on disk and returns only the
# HIV/AIDS disbursement rows, which are cached to FA_HIV_FILE; later runs read
# the small cache and never touch the big file.

if (!file.exists(FA_HIV_FILE)) {
  download_once(FA_URL, FA_FILE)
  con <- DBI::dbConnect(duckdb::duckdb())
  hiv_raw <- DBI::dbGetQuery(con, sprintf("
    SELECT \"Country Code\"          AS country_code,
           \"Country Name\"          AS country_name,
           \"Fiscal Year\"           AS fiscal_year,
           \"US Sector Name\"        AS us_sector_name,
           \"Managing Agency Name\"  AS managing_agency_name,
           \"Funding Account Name\"  AS funding_account_name,
           \"Activity Name\"         AS activity_name,
           \"Activity Description\"  AS activity_description,
           CAST(\"Current Dollar Amount\" AS DOUBLE) AS current_dollar_amount
    FROM read_csv('%s', all_varchar = true, header = true, quote = '\"',
                  escape = '\"', strict_mode = false, ignore_errors = true)
    WHERE \"Transaction Type Name\" = 'Disbursements'
      AND \"US Sector Name\" ILIKE '%%HIV%%'", FA_FILE))
  DBI::dbDisconnect(con, shutdown = TRUE)
  write_csv(hiv_raw, FA_HIV_FILE)
}

hiv <- read_csv(FA_HIV_FILE, show_col_types = FALSE,
                col_types = cols(.default = "c", current_dollar_amount = "d")) |>
  mutate(fiscal_year = suppressWarnings(as.integer(fiscal_year))) |>
  filter(!is.na(fiscal_year), fiscal_year <= max(YEARS))

# COVID-19 money is excluded two ways:
#  1. Activity name or funding account mentions COVID. Descriptions are NOT
#     searched: core PEPFAR awards (e.g. EpiC) mention COVID adaptations there.
#  2. HIV/AIDS disbursements from the Economic Support Fund in FY2021+. ESF
#     carried < $10M/yr of HIV money before 2021, then $1.5-2.1B in FY2021-22
#     (Global Fund contribution + PEPFAR partners) - the American Rescue Plan
#     COVID supplemental, which was routed through ESF.
covid_pattern <- "covid|coronavirus|sars-cov|ncov|american rescue plan|\\barpa\\b|cares act"
is_covid <- grepl(covid_pattern, paste(hiv$funding_account_name, hiv$activity_name),
                  ignore.case = TRUE) |
  (hiv$funding_account_name %in% "Economic Support Fund" & hiv$fiscal_year >= 2021)
message(sprintf("Dropping %d COVID-related HIV/AIDS rows ($%.1fM current).",
                sum(is_covid), sum(hiv$current_dollar_amount[is_covid], na.rm = TRUE) / 1e6))
hiv <- hiv[!is_covid, ]

# foreignassistance.gov country codes are ISO3; fall back to names for regional/odd codes
hiv_cy <- hiv |>
  mutate(iso3c = if_else(grepl("^[A-Z]{3}$", country_code), country_code, to_iso3(country_name))) |>
  filter(!is.na(iso3c)) |>
  group_by(iso3c, year = fiscal_year) |>
  summarise(us_hiv_cur = sum(current_dollar_amount, na.rm = TRUE), .groups = "drop")

pepfar <- hiv_cy |>
  filter(year >= PEPFAR_START, year %in% YEARS) |>
  transmute(iso3c, year, pepfar_cur = us_hiv_cur)

pre2004_dummy <- hiv_cy |>
  filter(year < PEPFAR_START, us_hiv_cur > 0) |>
  distinct(iso3c) |>
  mutate(us_hiv_pre2004 = 1L)

# =============================================================================
# 6. Country sample
# =============================================================================

countries <- read_excel(SAMPLE_FILE, sheet = "Countries") |>
  transmute(iso3c        = `ISO3 / WB code`,
            country      = `Country (WB name)`,
            study_group  = factor(Group, levels = c("COP-PEPFAR", "Other PEPFAR recipient", "Control")),
            grp_pepfar_any   = as.integer(pepfar_any),
            grp_pepfar_cop   = as.integer(pepfar_cop),
            grp_pepfar_other = as.integer(pepfar_other),
            grp_control      = as.integer(control))

# Kates-window-specific control flag -- see the note above "6. Country
# sample" for the full per-country threshold check. Identical to
# grp_control except VEN/COL/PAN (fail the Kates-window threshold) are 0.
KATES_CONTROL_EXCLUDE <- c("VEN", "COL", "PAN")
countries <- countries |>
  mutate(grp_control_kates = if_else(iso3c %in% KATES_CONTROL_EXCLUDE, 0L, grp_control),
         kates_threshold_note = if_else(
           iso3c %in% KATES_CONTROL_EXCLUDE,
           "Exceeds $1M AND $0.05/capita cumulative PEPFAR-linked disbursement 2004-2022; not a valid Kates-window control",
           ""))

sample_iso <- countries$iso3c
stopifnot(length(sample_iso) == N_EXPECTED, !anyDuplicated(sample_iso))

# Roughly half of HIV/AIDS disbursements each year are booked to "World" or a
# region (WLD, SSN, ...) rather than a country; those cannot enter the panel.
attributed <- hiv_cy |>
  filter(year >= PEPFAR_START) |>
  summarise(share = sum(us_hiv_cur[iso3c %in% sample_iso]) / sum(us_hiv_cur))
message(sprintf("Share of FY%d+ HIV/AIDS disbursements attributed to study countries: %.0f%%",
                PEPFAR_START, 100 * attributed$share))

# =============================================================================
# 7. OECD CRS: donor health ODA, excluding PEPFAR
# =============================================================================
# Health = DAC sectors 120 (Health) + 130 (Population/reproductive health).
# Non-PEPFAR = all-donor health minus US disbursements to purpose 13040
# (STD control incl. HIV/AIDS), which is where PEPFAR is reported in the CRS.

crs_get <- function(donor, sector) {
  key <- paste(donor, "", sector, "100", "_T", "_T", "D", "V", "_T", "", "", sep = ".")
  url <- paste0(CRS_BASE, "/", key,
                "?startPeriod=", min(YEARS), "&endPeriod=", max(YEARS),
                "&dimensionAtObservation=AllDimensions&format=csvfile")
  dest <- file.path(RAW_DIR, paste0("crs_", donor, "_", gsub("\\+", "-", sector), ".csv"))
  download_once(url, dest, mode = "wb")
  read_csv(dest, show_col_types = FALSE) |>
    filter(UNIT_MEASURE == "USD") |>
    # OBS_VALUE is scaled by UNIT_MULT (6 = millions)
    transmute(iso3c = recode(RECIPIENT, XKV = "XKX"),   # CRS codes Kosovo as XKV
              year = as.integer(TIME_PERIOD),
              value = as.numeric(OBS_VALUE) * 10^coalesce(as.numeric(UNIT_MULT), 0))
}

crs_health_all <- crs_get(CRS_DONOR, "120+130") |>
  group_by(iso3c, year) |> summarise(health_all = sum(value, na.rm = TRUE), .groups = "drop")

crs_us_hiv <- crs_get("USA", "13040") |>
  group_by(iso3c, year) |> summarise(us_hiv_crs = sum(value, na.rm = TRUE), .groups = "drop")

crs <- crs_health_all |>
  left_join(crs_us_hiv, by = c("iso3c", "year")) |>
  mutate(donor_health_nonpepfar_cur = pmax(health_all - coalesce(us_hiv_crs, 0), 0)) |>
  select(iso3c, year, donor_health_nonpepfar_cur)

# =============================================================================
# 8. GBD HIV prevalence (gap-fill for WDI)
# =============================================================================

if (file.exists(GBD_FILE)) {
  gbd <- read_csv(GBD_FILE, show_col_types = FALSE) |>
    rename_with(clean_names) |>
    filter(grepl("percent", metric_name, ignore.case = TRUE),
           grepl("15-49", age_name),
           grepl("both", sex_name, ignore.case = TRUE)) |>
    transmute(iso3c = to_iso3(location_name), year = as.integer(year),
              hiv_prev_gbd = val * 100) |>               # GBD "Percent" metric is a proportion
    filter(!is.na(iso3c))
} else {
  warning("GBD file not found; HIV prevalence will not be gap-filled: ", GBD_FILE)
  gbd <- tibble(iso3c = character(), year = integer(), hiv_prev_gbd = numeric())
}

# =============================================================================
# 9. Assemble panel
# =============================================================================

panel <- expand_grid(iso3c = sample_iso, year = YEARS) |>
  left_join(wdi,            by = c("iso3c", "year")) |>
  left_join(wpp,            by = c("iso3c", "year")) |>
  left_join(income_class,   by = c("iso3c", "year")) |>
  left_join(pepfar,         by = c("iso3c", "year")) |>
  left_join(crs,            by = c("iso3c", "year")) |>
  left_join(gbd,            by = c("iso3c", "year")) |>
  left_join(pre2004_dummy,  by = "iso3c") |>
  left_join(countries, by = "iso3c") |>
  mutate(
    us_hiv_pre2004 = coalesce(us_hiv_pre2004, 0L),

    # --- PEPFAR (0 in years with no disbursement; NA before program start)
    pepfar_cur = if_else(year >= PEPFAR_START, coalesce(pepfar_cur, 0), NA_real_),

    # --- HIV prevalence: WDI, gap-filled with GBD
    hiv_prev = coalesce(hiv_prev_wdi, hiv_prev_gbd),
    hiv_prev_source = case_when(!is.na(hiv_prev_wdi) ~ "WDI",
                                !is.na(hiv_prev_gbd) ~ "GBD",
                                TRUE ~ NA_character_),

    # --- Domestic health spending (government + private), PPP per capita
    dom_hexp_pc_ppp_cur = if_else(is.na(dom_hexp_gov_ppp) & is.na(dom_hexp_pvt_ppp), NA_real_,
                                  coalesce(dom_hexp_gov_ppp, 0) + coalesce(dom_hexp_pvt_ppp, 0)),

    # --- Constant 2022 USD
    gdp_pc               = deflate(gdp_pc_cur, year, defl),
    gdp                  = deflate(gdp_cur, year, defl),
    pepfar               = deflate(pepfar_cur, year, defl),
    donor_health_nonpepfar = deflate(donor_health_nonpepfar_cur, year, defl),
    dom_hexp_pc_ppp      = deflate(dom_hexp_pc_ppp_cur, year, defl),

    # --- Per capita (UN WPP population)
    pepfar_pc                 = pepfar / pop_total,
    donor_health_nonpepfar_pc = donor_health_nonpepfar / pop_total,
    pepfar_any                = as.integer(coalesce(pepfar, 0) > 0),

    # --- Primary school disengagement
    prim_disengage_girls = oos_prim_female / pop_prim_female,
    prim_disengage_boys  = oos_prim_male   / pop_prim_male,

    # --- Mortality: MMR only from 2000
    mort_maternal = if_else(year >= 2000, mort_maternal, NA_real_),

    # --- COVID period and health-burden interaction
    covid       = as.integer(year >= 2020),
    hiv_x_covid = hiv_prev * covid,

    # --- Analysis-window flags
    # FIX: in_kates was `year <= 2020`. Kates et al. (2026) covers 2004-2022,
    # so this now flags the paper's actual analysis window (and is therefore
    # 1 for every row, since the panel itself stops at 2022).
    in_crown = as.integer(year <= 2018),
    in_kates = as.integer(year <= 2022)
  )

# Ever-PEPFAR recipient and first PEPFAR year
panel <- panel |>
  group_by(iso3c) |>
  mutate(pepfar_ever = as.integer(any(pepfar_any == 1, na.rm = TRUE)),
         pepfar_first_year = if (any(pepfar_any == 1)) min(year[pepfar_any == 1]) else NA_integer_) |>
  ungroup()

# ---- Economic outcomes: levels, asinh, log(x+1) ----------------------------

econ_outcomes <- c("gdp_pc", "gdp", "gdp_pc_growth", "gdp_growth",
                   "emp_female", "emp_male",
                   "prim_disengage_girls", "prim_disengage_boys")

for (v in econ_outcomes) {
  panel[[paste0(v, "_ihs")]]   <- asinh(panel[[v]])
  panel[[paste0(v, "_log1p")]] <- safe_log1p(panel[[v]])   # NA where x <= -1 (e.g. negative growth)
}

# ---- Baseline covariates (Table 1), fixed at BASELINE_YR --------------------
# Uses the BASELINE_YR value; if missing, the latest non-missing value 2000-BASELINE_YR.

baseline_vars <- c(gdp_pc = "gdp_pc", pop_total = "pop_total", life_exp = "life_exp",
                   fertility = "fertility", urban_pct = "urban_pct",
                   sec_enroll_gross = "sec_enroll_gross", hiv_prev = "hiv_prev",
                   donor_health_nonpepfar_pc = "donor_health_nonpepfar_pc",
                   dom_hexp_pc_ppp = "dom_hexp_pc_ppp")

baseline <- panel |>
  filter(year >= 2000, year <= BASELINE_YR) |>
  arrange(iso3c, desc(year)) |>
  group_by(iso3c) |>
  summarise(across(all_of(unname(baseline_vars)), \(x) x[!is.na(x)][1]),
            bl_wb_income = wb_income[year == BASELINE_YR][1],
            .groups = "drop") |>
  rename_with(\(x) paste0("bl_", x), all_of(unname(baseline_vars)))

panel <- panel |>
  left_join(baseline, by = "iso3c") |>
  mutate(bl_us_hiv_pre2004 = us_hiv_pre2004)

# ---- Final column order -----------------------------------------------------

panel <- panel |>
  select(
    iso3c, country, year, study_group, starts_with("grp_"), kates_threshold_note,
    wb_income, covid, in_crown, in_kates,
    # treatment
    pepfar, pepfar_pc, pepfar_any, pepfar_ever, pepfar_first_year, us_hiv_pre2004,
    # economic outcomes
    all_of(econ_outcomes), ends_with("_ihs"), ends_with("_log1p"),
    # mortality outcomes
    mort_allcause, mort_child, mort_maternal,
    # time-varying covariates
    pop_total, life_exp, fertility, urban_pct, sec_enroll_gross,
    hiv_prev, hiv_prev_source, hiv_x_covid,
    donor_health_nonpepfar, donor_health_nonpepfar_pc, dom_hexp_pc_ppp,
    # baseline covariates
    starts_with("bl_"),
    # raw components kept for auditing
    gdp_pc_cur, gdp_cur, pepfar_cur, donor_health_nonpepfar_cur, dom_hexp_pc_ppp_cur,
    oos_prim_female, oos_prim_male, pop_prim_female, pop_prim_male,
    hiv_prev_wdi, hiv_prev_gbd
  ) |>
  arrange(iso3c, year)

# =============================================================================
# 10. Checks and export
# =============================================================================

stopifnot(!anyDuplicated(panel[c("iso3c", "year")]))
stopifnot(nrow(panel) == length(sample_iso) * length(YEARS))

message(sprintf("Panel: %d countries x %d years = %d rows",
                n_distinct(panel$iso3c), length(YEARS), nrow(panel)))

missing_report <- panel |>
  summarise(across(everything(), \(x) mean(is.na(x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "share_missing") |>
  arrange(desc(share_missing))
print(missing_report, n = 30)

codebook <- tribble(
  ~variable,                   ~description,                                              ~source,
  "study_group",               "COP-PEPFAR / Other PEPFAR recipient / Control",           "pepfar_country_list.xlsx",
  "grp_*",                     "Group dummies from the study country list",               "pepfar_country_list.xlsx",
  "grp_control_kates",         "Like grp_control, but VEN/COL/PAN set to 0 (fail the Kates-window $1M-and-$0.05pc control threshold, 2004-2022)", "derived",
  "kates_threshold_note",      "Why grp_control_kates differs from grp_control, where it does", "derived",
  "pepfar",                    "PEPFAR (US HIV/AIDS) disbursements, FY, excl. COVID, 2022 USD", "foreignassistance.gov",
  "pepfar_pc",                 "PEPFAR per capita, 2022 USD",                             "foreignassistance.gov / UN WPP",
  "us_hiv_pre2004",            "Received US HIV funding before FY2004 (dummy)",           "foreignassistance.gov",
  "gdp_pc",                    "GDP per capita, 2022 USD",                                "WDI NY.GDP.PCAP.CD (deflated)",
  "gdp",                       "GDP, 2022 USD",                                           "WDI NY.GDP.MKTP.CD (deflated)",
  "gdp_pc_growth",             "GDP per capita growth (annual %)",                        "WDI NY.GDP.PCAP.KD.ZG",
  "gdp_growth",                "GDP growth (annual %)",                                   "WDI NY.GDP.MKTP.KD.ZG",
  "emp_female",                "Employment-to-population ratio, female 15+ (%)",          "WDI SL.EMP.TOTL.SP.FE.ZS",
  "emp_male",                  "Employment-to-population ratio, male 15+ (%)",            "WDI SL.EMP.TOTL.SP.MA.ZS",
  "prim_disengage_girls",      "Out-of-school primary-age girls / primary-age girls",     "WDI SE.PRM.UNER.FE / SP.PRM.TOTL.FE.IN",
  "prim_disengage_boys",       "Out-of-school primary-age boys / primary-age boys",       "WDI SE.PRM.UNER.MA / SP.PRM.TOTL.MA.IN",
  "*_ihs / *_log1p",           "asinh(x) and log(x+1) of each economic outcome",          "derived",
  "mort_allcause",             "Crude death rate per 1,000",                              "WDI SP.DYN.CDRT.IN",
  "mort_child",                "Under-5 mortality per 1,000 live births",                 "WDI SH.DYN.MORT",
  "mort_maternal",             "Maternal mortality ratio per 100,000 live births (2000+)","WDI SH.STA.MMRT",
  "pop_total",                 "Total population (1 July)",                               "UN WPP 2024",
  "wb_income",                 "WB income group (L/LM/UM/H) for the data year",           "World Bank OGHIST",
  "life_exp",                  "Life expectancy at birth",                                "WDI SP.DYN.LE00.IN",
  "fertility",                 "Total fertility rate",                                    "WDI SP.DYN.TFRT.IN",
  "urban_pct",                 "Urban population (% of total)",                           "WDI SP.URB.TOTL.IN.ZS",
  "sec_enroll_gross",          "Secondary school enrollment (% gross)",                   "WDI SE.SEC.ENRR",
  "hiv_prev",                  "HIV prevalence, % ages 15-49 (WDI, GBD gap-fill)",        "WDI SH.DYN.AIDS.ZS / GBD",
  "hiv_x_covid",               "hiv_prev x covid (year >= 2020)",                         "derived",
  "donor_health_nonpepfar_pc", "Per capita donor health ODA excl. US HIV, 2022 USD",     "OECD CRS / UN WPP",
  "dom_hexp_pc_ppp",           "Domestic govt + private health exp per capita, PPP, 2022$","WDI SH.XPD.GHED.PP.CD + SH.XPD.PVTD.PP.CD",
  "bl_*",                      paste0("Baseline covariates at ", BASELINE_YR, " (latest 2000+ if missing)"), "derived"
)

write_csv(panel, file.path(OUT_DIR, "pepfar_panel_1990_2022.csv"), na = "")
saveRDS(panel,   file.path(OUT_DIR, "pepfar_panel_1990_2022.rds"))
write_csv(codebook, file.path(OUT_DIR, "codebook.csv"))
write_csv(countries, file.path(OUT_DIR, "sample_countries.csv"))
write_csv(missing_report, file.path(OUT_DIR, "missingness.csv"))

message("Done. Files written to ", normalizePath(OUT_DIR))

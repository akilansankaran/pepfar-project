# Step 3: share of PEPFAR dollars spent on maternal and child programs, by country and year
# Author: Rishad Vaghaiwalla
#
# WHY: cost-effectiveness for an outcome divides deaths averted by the dollars aimed at
# that outcome. Table 20 divides maternal deaths averted by ALL PEPFAR dollars, which can
# only give an upper bound on the cost per death. This script measures, for every PEPFAR
# country and fiscal year 2004-2022, the share of PEPFAR spending that went to
#   pmtct      prevention of mother-to-child transmission / pregnant and breastfeeding women
#   pediatric  pediatric care and treatment (children living with HIV)
#   ovc        orphans and vulnerable children
# life_years_setup.R multiplies these shares into the panel's PEPFAR disbursements
# (`pepfar`, constant 2022 USD) to get each outcome's matched dollars.
#
# SOURCES (all PEPFAR's own data):
#   1. PEPFAR Panorama Spotlight, "PEPFAR Operating Unit (OU) Budgets and Expenditures"
#      https://data.pepfar.gov/datasets (downloaded automatically below)
#        - expenditures FY2018-2022 by program, beneficiary and sub-beneficiary:
#          pmtct = beneficiary "Pregnant & Breastfeeding Women"; ovc = beneficiary "OVC";
#          pediatric = sub-beneficiary "Children", "Girls" or "Boys"
#        - budgets FY2004-2017 crosswalked to the new classification: only OVC is separable
#        - expenditure analysis FY2015-2017 (separate CSV in the same zip): program areas
#          "PMTCT" and "OVC"
#   2. PEPFAR Fiscal Year Operational Plans, FY2005-FY2011 (pepfar.gov, archived by the
#      UNT Libraries web archive), extracted into
#      data/intermediate/pepfar_operational_plans_FY2005_2011.csv (tracked in git):
#        - FY2010 and FY2011: approved funding by country and budget code
#          (MTCT = PMTCT; PDCS + PDTX = pediatric care and treatment; HKID = OVC)
#        - FY2005-FY2009: PEPFAR-wide totals by program area (no country detail published
#          as text; FY2009 is read from an image of Table 4)
#      Each row carries its source document, table and URL.
#
# OVERHEAD: management, strategic information and health-systems spending serve every
# program. Shares are computed "pro rata", i.e. as a share of program spending excluding
# that overhead, which is the same as giving each program its proportional slice of the
# overhead. (FY2006-FY2008 PEPFAR already attributed overhead to program areas this way.)
# Overhead: PM and ASP (expenditure reporting); PM, HSS, SI, SURV (expenditure analysis);
# M&O, HVSI, OHSS (operational plans).
#
# FILLING GAPS: each country gets its own share where a source covers it. In other years,
# its nearest observed share is moved along the PEPFAR-wide trend:
#   share_c,t = share_c,t* x G_t / G_t*   (t* = nearest year with country data)
# where G_t is the PEPFAR-wide share (published totals, or the dollar-weighted average of
# the country data; linear interpolation between source years). Countries with no own data
# get G_t. Every value records where it came from (src_* columns).
#
# DEFINITION BREAK: from FY2018 PEPFAR records PMTCT as treatment and testing spending
# "targeted" at pregnant and breastfeeding women; ART for pregnant women in general
# clinics is no longer separable, so the PMTCT share falls (roughly 6-7% to 2%).
# pmtct_after_2017 = "er" uses the new definition as reported; "hold_ea" carries each
# country's FY2015-2017 PMTCT share forward instead (sensitivity).
#
# Run from the repository root:  Rscript data_pull/build_pepfar_outcome_shares.R
# Writes data/intermediate/pepfar_outcome_shares.csv and, in
# tables/cost_effectiveness/by_outcome/: csv/pepfar_spending_shares_by_year.csv and
# figures/pepfar_spending_shares.{pdf,png}

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(countrycode); library(ggplot2)
})

# =====================================================
# 0. Settings
# =====================================================
spend_dir <- Sys.getenv("PEPFAR_SPENDING_DIR", "data/raw/pepfar_spending")
ou_url    <- "https://data.pepfar.gov/media/fin_mngt/PEPFAR_OU_Budgets_and_Expenditures.zip"
ou_zip    <- file.path(spend_dir, "PEPFAR_OU_Budgets_and_Expenditures.zip")
er_file   <- file.path(spend_dir, "PEPFAR_OU_Budgets_and_Program_Expenditures.txt")
ea_file   <- file.path(spend_dir, "PEPFAR_OU_Expenditures_FY_2015_2017.csv")
op_file   <- "data/intermediate/pepfar_operational_plans_FY2005_2011.csv"   # extracted tables, tracked in git
out_file  <- "data/intermediate/pepfar_outcome_shares.csv"
ce_dir    <- "tables/cost_effectiveness/by_outcome"
years     <- 2004:2022

# Which outcome PMTCT dollars count toward. PMTCT serves mothers (ART during pregnancy and
# breastfeeding) and infants (preventing infection), so the default counts it in both.
# Each outcome's cost-effectiveness is computed separately, so this double counts nothing
# within an outcome; just don't add maternal and child dollars together.
#   "both" (default), "maternal", "child", or "split" (half each)
pmtct_to <- "both"
pmtct_after_2017 <- "er"   # "er" (reported definition) or "hold_ea" (see DEFINITION BREAK)

# =====================================================
# 1. Download PEPFAR's OU budgets and expenditures (same pattern as the WPP downloads)
# =====================================================
dir.create(spend_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(er_file) || !file.exists(ea_file)) {
  if (!file.exists(ou_zip)) {
    message("Downloading ", basename(ou_zip), " from data.pepfar.gov ...")
    options(timeout = max(600, getOption("timeout")))
    ok <- tryCatch({ download.file(ou_url, ou_zip, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) stop("Could not download ", ou_url, ". Download it in a browser from ",
                  "https://data.pepfar.gov/datasets (Financial Management Datasets) and save it to ", spend_dir)
  }
  unzip(ou_zip, exdir = spend_dir)
}
if (!file.exists(op_file)) stop(op_file, " is missing (it comes with these scripts; see its source and url columns).")

iso_of <- function(x) {
  fix <- c("CotedIvoire" = "CIV", "Cote d' Ivoire" = "CIV", "DRC" = "COD",
           "Democratic Republic of the" = "COD", "Burma" = "MMR", "SouthSudan" = "SSD")
  out <- unname(fix[x])
  need <- is.na(out)
  out[need] <- suppressWarnings(countrycode(gsub("([a-z])([A-Z])", "\\1 \\2", x[need]),
                                            "country.name", "iso3c"))
  out
}

# =====================================================
# 2. Expenditure reporting (FY2018-2022) and budgets (FY2004-2017, OVC only)
# =====================================================
er <- read_tsv(er_file, col_types = cols(.default = col_character()), progress = FALSE) |>
  rename_with(~ gsub("^﻿", "", .x)) |>
  filter(`Operating Unit` != "Total", !Program %in% c("Total", "Applied Pipeline"), !is.na(Program),
         !is.na(ISO3), ISO3 != "") |>
  pivot_longer(matches("^\\d{4} (BUDGET|EXPENDITURE)$"), names_to = c("year", "kind"),
               names_sep = " ", values_to = "usd") |>
  mutate(year = as.integer(year), usd = parse_number(usd)) |>
  filter(!is.na(usd), usd != 0) |>
  mutate(bucket = case_when(
    Program %in% c("PM", "ASP") ~ "overhead",
    Beneficiary %in% "Pregnant & Breastfeeding Women" | `Targeted Beneficiary` %in% "Pregnant & Breastfeeding Women" ~ "pmtct",
    Beneficiary %in% "OVC" | `Targeted Beneficiary` %in% "OVC" ~ "ovc",
    `Sub Beneficiary` %in% c("Children", "Girls", "Boys") | `Targeted Beneficiary` %in% "Children" ~ "pediatric",
    TRUE ~ "other"))

shares_from <- function(d, iso, src) {
  d |>
    group_by(iso3c = {{ iso }}, year) |>
    summarise(denom = sum(usd[bucket != "overhead"]), pmtct = sum(usd[bucket == "pmtct"]),
              pediatric = sum(usd[bucket == "pediatric"]), ovc = sum(usd[bucket == "ovc"]),
              .groups = "drop") |>
    filter(denom > 0) |>
    mutate(source = src)
}
er_obs  <- er |> filter(kind == "EXPENDITURE", year %in% 2018:2022) |> shares_from(ISO3, "PEPFAR expenditure reporting")
# FY2004 budgets have no program detail, so OVC budget shares start in FY2005
bud_obs <- er |> filter(kind == "BUDGET", year %in% 2005:2017) |> shares_from(ISO3, "PEPFAR OU budgets (crosswalked)")

# =====================================================
# 3. Expenditure analysis (FY2015-2017)
# =====================================================
ea_txt <- readLines(ea_file, warn = FALSE, encoding = "UTF-8")   # the file uses CR line endings
ea <- read_csv(I(ea_txt), col_types = cols(.default = col_character()), progress = FALSE) |>
  rename_with(~ gsub("^﻿", "", .x)) |>
  pivot_longer(starts_with("Fiscal Year"), names_to = "year", values_to = "usd") |>
  mutate(year = as.integer(sub("Fiscal Year ", "", year)), usd = parse_number(usd),
         iso3c = iso_of(`Benefiting Country`),
         bucket = case_when(`Expenditure Analysis Program Area` == "PMTCT" ~ "pmtct",
                            `Expenditure Analysis Program Area` == "OVC" ~ "ovc",
                            `Expenditure Analysis Program Area` %in% c("PM", "HSS", "SI", "SURV") ~ "overhead",
                            TRUE ~ "other")) |>
  filter(!is.na(usd), usd != 0, !is.na(iso3c))
ea_obs <- shares_from(ea, iso3c, "PEPFAR expenditure analysis") |> mutate(pediatric = NA_real_)

# =====================================================
# 4. Operational plans (FY2005-2011)
# =====================================================
op <- read_csv(op_file, show_col_types = FALSE) |>
  mutate(denom = total - coalesce(overhead, 0), ped = ped_care + ped_treatment)
op_obs <- op |> filter(level == "country") |>
  transmute(iso3c = iso_of(operating_unit), year = fiscal_year, denom, pmtct, pediatric = ped, ovc,
            source = "PEPFAR operational plan") |>
  filter(!is.na(iso3c))
op_glob <- op |> filter(level == "global") |>
  transmute(year = fiscal_year, pmtct = pmtct / denom, pediatric = ped / denom, ovc = ovc / denom)

# =====================================================
# 5. Observed country-year shares, one source per bucket
# =====================================================
obs_long <- bind_rows(
  op_obs |> transmute(iso3c, year, denom, bucket = "pmtct", usd = pmtct, source),
  op_obs |> transmute(iso3c, year, denom, bucket = "pediatric", usd = pediatric, source),
  ea_obs |> transmute(iso3c, year, denom, bucket = "pmtct", usd = pmtct, source),
  er_obs |> transmute(iso3c, year, denom, bucket = "pmtct", usd = pmtct, source),
  er_obs |> transmute(iso3c, year, denom, bucket = "pediatric", usd = pediatric, source),
  bud_obs |> transmute(iso3c, year, denom, bucket = "ovc", usd = ovc, source),
  er_obs |> transmute(iso3c, year, denom, bucket = "ovc", usd = ovc, source)) |>
  filter(!is.na(usd), denom > 0) |>
  mutate(share = pmin(usd / denom, 1))
if (pmtct_after_2017 == "hold_ea") obs_long <- filter(obs_long, !(bucket == "pmtct" & year >= 2018))

# PEPFAR-wide series G_t: published totals where they exist, else the dollar-weighted
# average of the country data, else linear interpolation (flat at the ends)
agg <- obs_long |> group_by(bucket, year) |> summarise(g = sum(usd) / sum(denom), .groups = "drop")
glob <- expand_grid(bucket = c("pmtct", "pediatric", "ovc"), year = years) |>
  left_join(op_glob |> pivot_longer(-year, names_to = "bucket", values_to = "g_pub"), by = c("bucket", "year")) |>
  left_join(agg, by = c("bucket", "year")) |>
  mutate(g = coalesce(g_pub, g),
         g_source = case_when(!is.na(g_pub) ~ "operational plan total", !is.na(g) ~ "country data, dollar-weighted",
                              TRUE ~ "interpolated")) |>
  group_by(bucket) |>
  mutate(g = approx(year[!is.na(g)], g[!is.na(g)], xout = year, rule = 2)$y) |>
  ungroup() |> select(bucket, year, g, g_source)
if (pmtct_after_2017 == "hold_ea") glob <- glob |>
  group_by(bucket) |> mutate(g = ifelse(bucket == "pmtct" & year >= 2018, g[year == 2017], g)) |> ungroup()

# =====================================================
# 6. Every PEPFAR-funded country-year in the panel gets a share
# =====================================================
panel <- read_csv("data/intermediate/pepfar_panel_1990_2022_FIXED.csv",
                  col_select = c(iso3c, year, pepfar), show_col_types = FALSE) |>
  filter(year %in% years, coalesce(pepfar, 0) != 0)

fill_bucket <- function(b) {
  o <- filter(obs_long, bucket == b) |> select(iso3c, year, share, source)
  g <- filter(glob, bucket == b) |> select(year, g)
  panel |> select(iso3c, year) |>
    left_join(o, by = c("iso3c", "year")) |>
    left_join(g, by = "year") |>
    group_by(iso3c) |>
    group_modify(function(d, key) {
      own <- filter(o, iso3c == key$iso3c)
      if (!nrow(own)) return(mutate(d, share = g, source = "PEPFAR-wide share (no country data)"))
      gy <- setNames(g$g, g$year)
      for (i in which(is.na(d$share))) {
        j <- which.min(abs(own$year - d$year[i]) + 0.1 * (own$year > d$year[i]))  # ties: earlier year
        d$share[i]  <- min(1, own$share[j] * gy[[as.character(d$year[i])]] / gy[[as.character(own$year[j])]])
        d$source[i] <- paste0("own ", own$year[j], " share x PEPFAR-wide trend")
      }
      d
    }) |>
    ungroup() |>
    transmute(iso3c, year, !!paste0("share_", b) := share, !!paste0("src_", b) := source)
}
out <- Reduce(function(a, b) full_join(a, b, by = c("iso3c", "year")),
              lapply(c("pmtct", "pediatric", "ovc"), fill_bucket))

stopifnot(pmtct_to %in% c("both", "maternal", "child", "split"))
w_mat   <- switch(pmtct_to, both = 1, maternal = 1, split = 0.5, child = 0)
w_child <- switch(pmtct_to, both = 1, child = 1, split = 0.5, maternal = 0)
out <- out |>
  mutate(share_maternal = share_pmtct * w_mat,
         share_child    = pmin(1, share_pediatric + share_ovc + share_pmtct * w_child)) |>
  relocate(iso3c, year, share_maternal, share_child)
write_csv(out, out_file)

# =====================================================
# 7. Summary by year and a figure
# =====================================================
by_year <- out |>
  left_join(panel, by = c("iso3c", "year")) |>
  group_by(year) |>
  summarise(pepfar_usd_m = sum(pepfar) / 1e6,
            across(c(share_pmtct, share_pediatric, share_ovc, share_maternal, share_child),
                   ~ sum(.x * pepfar) / sum(pepfar)),
            pct_pmtct_own = 100 * mean(src_pmtct %in% unique(obs_long$source)),
            pct_ped_own   = 100 * mean(src_pediatric %in% unique(obs_long$source)), .groups = "drop") |>
  left_join(glob |> select(-g) |> pivot_wider(names_from = bucket, values_from = g_source, names_prefix = "g_src_"),
            by = "year")
dir.create(file.path(ce_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ce_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
write_csv(by_year, file.path(ce_dir, "csv", "pepfar_spending_shares_by_year.csv"))

pal <- c("#2a78d6", "#eb6834", "#1baf7a")
fig <- glob |>
  mutate(bucket = recode(bucket, pmtct = "PMTCT (mothers and infants)", pediatric = "Pediatric care and treatment",
                         ovc = "Orphans and vulnerable children")) |>
  ggplot(aes(year, 100 * g, color = bucket)) +
  geom_vline(xintercept = 2017.5, color = "#52514e", linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.9) + geom_point(aes(shape = g_source), size = 1.8) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c("operational plan total" = 15, "country data, dollar-weighted" = 16, "interpolated" = 1)) +
  labs(title = "Share of PEPFAR program spending by target group, PEPFAR-wide",
       subtitle = "Overhead allocated pro rata. Dashed: FY2018 switch to expenditure reporting (narrower PMTCT definition).",
       x = "Fiscal year", y = "% of program spending", color = NULL, shape = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", legend.box = "vertical", legend.justification = "left",
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(ce_dir, "figures", "pepfar_spending_shares.pdf"), fig, width = 8, height = 4.8)
ggsave(file.path(ce_dir, "figures", "pepfar_spending_shares.png"), fig, width = 8, height = 4.8, dpi = 200, bg = "white")

cat("\nPMTCT dollars counted toward:", pmtct_to, "| PMTCT after FY2017:", pmtct_after_2017, "\n")
cat("\nPEPFAR-wide shares (G_t) and their source:\n")
print(glob |> mutate(g = round(100 * g, 1)) |> pivot_wider(names_from = bucket, values_from = c(g, g_source)), n = Inf, width = Inf)
cat("\nDollar-weighted shares applied to the panel's PEPFAR disbursements:\n")
print(by_year |> transmute(year, pepfar_usd_m = round(pepfar_usd_m),
                           across(c(share_pmtct, share_pediatric, share_ovc, share_maternal, share_child), ~ round(100 * .x, 1)),
                           pct_pmtct_own = round(pct_pmtct_own), pct_ped_own = round(pct_ped_own)), n = Inf, width = Inf)
cat("\nWrote", out_file, "\n")

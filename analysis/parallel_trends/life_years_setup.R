# Shared helpers for the life-years and cost-effectiveness scripts (maternal, child, all-cause)
# Author: Rishad Vaghaiwalla
#
# This is the machinery of life_years_maternal.R, lifted out so the child and all-cause
# versions are built exactly the same way: same trend-break model, same counterfactual,
# same WPP 2024 inputs, same discounting, same clustered-SE totals. Only the outcome, its
# denominator (births or population) and whose life expectancy a death averted is worth
# change, and those are set in `outcome_cfg` below.
#
# Check: with y = "mort_maternal", run_scenario() reproduces
# tables/cost_effectiveness/csv/maternal_life_years_scenarios.csv (life_years_maternal.R).
#
# Sourced by:
#   analysis/parallel_trends/life_years_mortality.R          (Step 2: Table 20 for each outcome)
#   analysis/parallel_trends/counterfactual_plots_mortality.R (Step 2: counterfactual plots)
#   analysis/parallel_trends/cost_benefit_income_groups.R    (Step 5: CE and B/C by income group)
#
# Model (Burgess and Pande 2005, eq. 3; see trend_break_mortality.R):
#   Y_it = a_i + g_t + γ1 (X_i x [t - 1990]) + γ2 (X_i x [t - B] x Post_t(B))
#                    + γ4 (X_i x Post_t(B)) + e_it
# Counterfactual "trend": γ2 = γ4 = 0, so for t >= B the gap is
#   Y_cf - Y_obs = -X_i (γ4 + γ2 (t - B))           [deaths per `scale` of the denominator]
# "flat": the pre-break trend also stops at B:  -X_i (γ4 + (γ1 + γ2)(t - B))

source("analysis/econ_setup.R")

# =====================================================
# 0. Settings shared by every outcome
# =====================================================
start_year <- 1990
last_year  <- 2022
r_main     <- 0.03          # discount rate
spec_main  <- "Adjusted"    # main specification
version_main <- "HIV burden"  # X_i = HIV prevalence 2003, as in Table 20

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
ce_dir  <- "tables/cost_effectiveness/by_outcome"
sweep_best_file <- "tables/kates_trend_breaks/sweep/csv/best_break_year.csv"  # from Step 4
shares_file     <- "data/intermediate/pepfar_outcome_shares.csv"              # from Step 3

# Outcome settings
#   scale      the rate is per `scale` units of the denominator
#   denom      "births" (WPP Births) or "pop" (WPP mid-year population)
#   sex, age   remaining life expectancy at this exact age (WPP abridged life table:
#              ages 0, 1, 5, 10, ..., 100) is what one death averted is worth
#   ages_alt   the two "age at death" robustness rows
#   break_year main break year; NA = the best-fitting year from Step 4
#              (tables/kates_trend_breaks/sweep/csv/best_break_year.csv, HIV burden, adjusted)
#   share_col  column of data/intermediate/pepfar_outcome_shares.csv (Step 3) with the share
#              of PEPFAR dollars matched to this outcome; NA = all PEPFAR dollars
outcome_cfg <- list(
  mort_maternal = list(label = "Maternal", file = "maternal", long = "maternal deaths",
                       rate_tex = "maternal deaths per 100,000 live births",
                       scale = 1e5, denom = "births", sex = "Female", age = 25,
                       ages_alt = c(20, 30), le_note = "remaining female life expectancy",
                       break_year = 2008, share_col = "share_maternal"),
  mort_child    = list(label = "Under-5", file = "child", long = "under-5 deaths",
                       rate_tex = "under-5 deaths per 1,000 live births",
                       scale = 1e3, denom = "births", sex = "Total", age = 0,
                       ages_alt = c(1, 5), le_note = "life expectancy at birth (both sexes)",
                       break_year = NA, share_col = "share_child"),
  mort_allcause = list(label = "All-cause", file = "allcause", long = "deaths (all causes)",
                       rate_tex = "deaths per 1,000 people (crude death rate)",
                       scale = 1e3, denom = "pop", sex = "Total", age = 35,
                       ages_alt = c(25, 45), le_note = "remaining life expectancy (both sexes)",
                       break_year = NA, share_col = NA)
)

main_break <- function(y) {
  B <- outcome_cfg[[y]]$break_year
  if (!is.na(B)) return(B)
  if (!file.exists(sweep_best_file))
    stop("No break year set for ", y, ". Run analysis/parallel_trends/trend_break_sweep_mortality.R",
         " first (it writes ", sweep_best_file, "), or set break_year in outcome_cfg.")
  b <- read_csv(sweep_best_file, show_col_types = FALSE) |>
    filter(outcome == y, version == version_main, spec == spec_main)
  stopifnot(nrow(b) == 1)
  b$best_break
}

# =====================================================
# 1. Trend-break model (same as trend_break_maternal_2008.R, any outcome and version)
# =====================================================
versions <- list(
  "All PEPFAR" = list(data = samples[["All PEPFAR"]], xvar = "treat"),
  "COP"        = list(data = samples[["COP"]],        xvar = "treat"),
  "HIV burden" = list(data = samples[["All PEPFAR"]], xvar = "bl_hiv_prev")
)

add_terms <- function(d, xvar, B) {
  d |>
    filter(!is.na(.data[[xvar]])) |>
    mutate(x         = .data[[xvar]],
           post_b    = as.integer(year >= B),
           tr_base   = x * (year - start_year),   # γ1
           tr_change = x * (year - B) * post_b,   # γ2
           jump      = x * post_b)                # γ4
}
fe_part <- function(spec, xvar) {
  if (spec == "Unadjusted") return("iso3c + year")
  ctrl <- if (xvar == "bl_hiv_prev") setdiff(bl_cont, "bl_hiv_prev") else bl_cont
  paste0("iso3c + bl_wb_income^year + year[", paste(ctrl, collapse = ", "), "]")
}
fit_cache <- list()
fit_model <- function(y, B, spec = spec_main, version = version_main) {
  key <- paste(y, B, spec, version)
  if (is.null(fit_cache[[key]])) {
    vv <- versions[[version]]
    d  <- add_terms(vv$data, vv$xvar, B)
    m  <- feols(as.formula(paste0(y, " ~ tr_base + tr_change + jump | ", fe_part(spec, vv$xvar))),
                data = d, cluster = ~iso3c)
    fit_cache[[key]] <<- list(m = m, used = d[obs(m), ])
  }
  fit_cache[[key]]
}

# =====================================================
# 2. WPP 2024: births, population, life expectancy (same files as life_years_maternal.R)
# =====================================================
wpp_base <- "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/"
di_file  <- file.path(raw_dir, "WPP2024_Demographic_Indicators_Medium.csv.gz")
lt_file  <- file.path(raw_dir, "WPP2024_Life_Table_Abridged_Medium_1950-2023.csv.gz")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

get_file <- function(file) {
  if (!file.exists(file)) {
    message("Downloading ", basename(file), " ...")
    options(timeout = max(3600, getOption("timeout")))
    ok <- tryCatch({ download.file(paste0(wpp_base, basename(file)), file, mode = "wb", quiet = TRUE)
      TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok && file.exists(file)) file.remove(file)
  }
  file.exists(file)
}
for (f in c(di_file, lt_file)) if (!get_file(f))
  stop("Could not download ", basename(f), ". Download it from the WPP 2024 CSV page and save it to ", raw_dir)

di <- read_csv(di_file, col_select = c(ISO3_code, Time, Variant, Births, TPopulation1July, LEx),
               show_col_types = FALSE) |>
  filter(!is.na(ISO3_code), Variant == "Medium", Time >= 2000, Time <= last_year) |>
  transmute(iso3c = ISO3_code, year = as.integer(Time),
            births = Births * 1000,                 # WPP reports thousands
            pop    = TPopulation1July * 1000,
            le0    = LEx) |>
  distinct(iso3c, year, .keep_all = TRUE)

lt <- read_csv(lt_file, col_select = c(ISO3_code, Time, Sex, AgeGrpStart, ex), show_col_types = FALSE) |>
  filter(!is.na(ISO3_code), Sex %in% c("Female", "Total"), Time >= 2000, Time <= last_year) |>
  transmute(iso3c = ISO3_code, year = as.integer(Time), sex = Sex, age = AgeGrpStart, ex) |>
  distinct(iso3c, year, sex, age, .keep_all = TRUE)
ex_source <- "WPP 2024 abridged life table"

make_ex <- function(age_x, sex_x) {
  if (!age_x %in% lt$age) stop("age must be one of: ", paste(sort(unique(lt$age)), collapse = ", "))
  lt |> filter(age == age_x, sex == sex_x) |> select(iso3c, year, ex)
}

# =====================================================
# 3. Matched PEPFAR dollars (Step 3). Without the shares file, all PEPFAR dollars are used,
#    exactly as in Table 20.
# =====================================================
shares <- if (file.exists(shares_file)) read_csv(shares_file, show_col_types = FALSE) else NULL
if (is.null(shares)) message("Note: ", shares_file, " not found; maternal and child costs use ALL ",
                             "PEPFAR dollars (run data_pull/build_pepfar_outcome_shares.R first).")
dollars_note <- function(y) {
  sc <- outcome_cfg[[y]]$share_col
  if (is.na(sc)) "all PEPFAR disbursements"
  else if (is.null(shares)) "all PEPFAR disbursements (outcome shares not yet available, so an upper bound)"
  else paste0("PEPFAR disbursements times the ", outcome_cfg[[y]]$label, " share of program spending")
}

# =====================================================
# 4. Counterfactual rows and totals
# =====================================================
# Country-year rows with every input, t >= B; k = years since the break.
# scope "sample": countries in the model's estimation sample; "all": also applies the
# coefficients to countries the adjusted model drops for missing baseline covariates.
cy_rows <- function(y, B, spec = spec_main, version = version_main, age_x = NULL, r = r_main,
                    ex_adj = 1, scope = "sample") {
  cfg <- outcome_cfg[[y]]
  if (is.null(age_x)) age_x <- cfg$age
  vv  <- versions[[version]]
  src <- if (scope == "all") filter(add_terms(vv$data, vv$xvar, B), !is.na(.data[[y]])) else
    fit_model(y, B, spec, version)$used
  out <- src |>
    filter(year >= B) |>
    select(iso3c, country, year, treat, x, all_of(y), gdp_pc, pepfar, wb_income) |>
    left_join(select(di, iso3c, year, denom = all_of(cfg$denom), le0), by = c("iso3c", "year")) |>
    left_join(make_ex(age_x, cfg$sex), by = c("iso3c", "year")) |>
    mutate(k      = year - B,
           ex_use = ex * ex_adj,
           ann    = if (r > 0) (1 - (1 + r)^(-ex_use)) / r else ex_use,  # discounted years per death
           disc   = (1 + r)^(-k),                                       # back to the break year
           pepfar = coalesce(pepfar, 0))
  # matched dollars: share x PEPFAR (share = 1 for all-cause, or when shares are unavailable)
  if (!is.na(cfg$share_col) && !is.null(shares)) {
    out <- out |>
      left_join(select(shares, iso3c, year, share = all_of(cfg$share_col)), by = c("iso3c", "year")) |>
      mutate(pepfar_match = pepfar * coalesce(share, 0))
  } else {
    out <- mutate(out, pepfar_match = pepfar)
  }
  out
}

# Weights on (γ1, γ2, γ4) that turn the coefficients into deaths averted in each row
row_w <- function(d, cf, jump, scale) {
  base <- d$x * d$denom / scale
  cbind(tr_base   = if (cf == "flat") -base * d$k else 0,
        tr_change = -base * d$k,
        jump      = if (jump) -base else 0)
}

# Total of f x (deaths averted) over the rows, with its clustered SE (w' V w)
lin_total <- function(W, f, b, V) {
  w <- colSums(W * f)
  nm <- names(w)
  c(est = sum(w * b[nm]), se = sqrt(drop(t(w) %*% V[nm, nm] %*% w)))
}

# Totals for one set of rows. `value_per_ly` is the dollar value of one discounted life-year
# in each row (GDP per capita in Table 20; a threshold or a VSLY in Step 5).
totals <- function(dd, W, b, V, value_per_ly = dd$gdp_pc) {
  deaths <- lin_total(W, 1, b, V)
  ly_u   <- lin_total(W, dd$ex_use, b, V)
  ly_d   <- lin_total(W, dd$ann * dd$disc, b, V)
  val    <- lin_total(W, dd$ann * dd$disc * value_per_ly, b, V)
  tibble(n_countries = n_distinct(dd$iso3c), country_years = nrow(dd),
         deaths_averted = deaths[["est"]], deaths_se = deaths[["se"]],
         ly_undisc = ly_u[["est"]], ly_undisc_se = ly_u[["se"]],
         ly_disc   = ly_d[["est"]], ly_disc_se   = ly_d[["se"]],
         value_usd = val[["est"]],  value_se     = val[["se"]],
         pepfar_disc_usd       = sum(dd$pepfar * dd$disc),
         pepfar_match_disc_usd = sum(dd$pepfar_match * dd$disc)) |>
    mutate(deaths_lo = deaths_averted - 1.96 * deaths_se, deaths_hi = deaths_averted + 1.96 * deaths_se,
           ly_lo     = ly_disc - 1.96 * ly_disc_se,       ly_hi     = ly_disc + 1.96 * ly_disc_se,
           # all PEPFAR dollars (as in Table 20)
           cost_per_ly    = ifelse(ly_disc > 0, pepfar_disc_usd / ly_disc, NA_real_),
           cost_per_ly_lo = ifelse(ly_hi > 0, pepfar_disc_usd / ly_hi, NA_real_),  # more LY -> cheaper
           cost_per_ly_hi = ifelse(ly_lo > 0, pepfar_disc_usd / ly_lo, NA_real_),
           benefit_cost   = value_usd / pepfar_disc_usd,
           # matched dollars (Step 3); cost-effectiveness = dollars per death averted
           cost_per_death_match = ifelse(deaths_averted > 0, pepfar_match_disc_usd / deaths_averted, NA_real_),
           cost_per_ly_match    = ifelse(ly_disc > 0, pepfar_match_disc_usd / ly_disc, NA_real_),
           benefit_cost_match   = value_usd / pepfar_match_disc_usd)
}

run_scenario <- function(y, label, B = main_break(y), spec = spec_main, version = version_main,
                         age_x = NULL, r = r_main, ex_adj = 1, cf = "trend", jump = TRUE,
                         scope = "sample") {
  cfg <- outcome_cfg[[y]]
  if (is.null(age_x)) age_x <- cfg$age
  fm <- fit_model(y, B, spec, version); b <- coef(fm$m); V <- vcov(fm$m)
  all_rows <- cy_rows(y, B, spec, version, age_x, r, ex_adj, scope)
  ok <- complete.cases(all_rows[, c("denom", "ex", "gdp_pc")])
  d  <- all_rows[ok, ]
  bind_rows(lapply(c("All countries", "PEPFAR recipients"), function(g) {
    dd <- if (g == "PEPFAR recipients") filter(d, treat == 1) else d
    tibble(outcome = y, scenario = label, group = g, break_year = B, spec = spec, version = version,
           cf = cf, jump = jump, scope = scope, age_death = age_x, r = r, ex_adj = ex_adj,
           dropped_missing = sum(!ok),
           obs_deaths = sum(dd[[y]] * dd$denom / cfg$scale)) |>
      bind_cols(totals(dd, row_w(dd, cf, jump, cfg$scale), b, V))
  }))
}

# Main-scenario rows with point estimates per country-year (for by-country, by-year and
# income-group breakdowns). Returns the rows plus b and V so totals can carry SEs.
main_rows <- function(y, B = main_break(y), spec = spec_main, version = version_main) {
  cfg <- outcome_cfg[[y]]
  fm <- fit_model(y, B, spec, version); b <- coef(fm$m); V <- vcov(fm$m)
  d  <- cy_rows(y, B, spec, version) |> filter(complete.cases(denom, ex, gdp_pc))
  W  <- row_w(d, "trend", TRUE, cfg$scale)
  d  <- d |> mutate(deaths_obs     = .data[[y]] * denom / cfg$scale,
                    deaths_averted = drop(W %*% b[colnames(W)]),
                    deaths_cf      = deaths_obs + deaths_averted,
                    ly_disc        = deaths_averted * ann * disc,
                    value_usd      = ly_disc * gdp_pc)
  list(rows = d, W = W, b = b, V = V, B = B)
}

# Observed vs. counterfactual deaths per year, treated countries, with a 95% band
by_year_cf <- function(mr) {
  d <- mutate(mr$rows, .i = row_number()) |> filter(treat == 1)
  bind_rows(lapply(sort(unique(d$year)), function(t) {
    dd  <- filter(d, year == t)
    est <- lin_total(mr$W[dd$.i, , drop = FALSE], 1, mr$b, mr$V)
    tibble(year = t, deaths_obs = sum(dd$deaths_obs),
           deaths_averted = est[["est"]], deaths_averted_se = est[["se"]])
  })) |>
    mutate(deaths_cf = deaths_obs + deaths_averted,
           cf_lo = deaths_cf - 1.96 * deaths_averted_se,
           cf_hi = deaths_cf + 1.96 * deaths_averted_se)
}

# =====================================================
# 5. Output helpers (same formats as life_years_maternal.R)
# =====================================================
save_both <- function(p, dir, name, w, h) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
fmt0   <- function(v) formatC(round(v), format = "d", big.mark = ",")
fmt_ci <- function(lo, hi) sprintf("[%s, %s]", fmt0(lo), fmt0(hi))
fmt_bn <- function(v) sprintf("%.2f", v / 1e9)
big_lab <- function(v) formatC(v, format = "d", big.mark = ",")

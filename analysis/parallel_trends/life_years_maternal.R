# Life-years saved: maternal mortality counterfactual from the trend-break model
# Author: Rishad Vaghaiwalla
#
# Takes the maternal mortality trend-break model (HIV burden version, adjusted, break in
# B = 2008; see trend_break_maternal_2008.R) and asks how many maternal deaths and life-years
# the post-break change is associated with.
#
# Counterfactual ("trend"): no break. Each country stays on its pre-break path, i.e. the
# fitted model with γ2 = γ4 = 0. For t >= B the gap between counterfactual and observed
# maternal mortality is
#     MMR_cf - MMR_obs = -X_i (γ4 + γ2 (t - B))          [deaths per 100,000 live births]
# Alternative counterfactual ("flat"): the pre-break trend stops at B instead of continuing,
#     MMR_cf - MMR_obs = -X_i (γ4 + (γ1 + γ2)(t - B))
# The "flat" version matters because γ1 > 0 here: the "trend" counterfactual assumes
# maternal mortality in high-HIV countries would have kept rising after B.
#
# Then, per country-year t >= B:
#   deaths averted   = gap x births / 100,000                       (WPP 2024 births)
#   life-years       = deaths averted x remaining female life expectancy at age_death
#                      (WPP 2024 life table), discounted at r within each life and
#                      back to B across years
#   value            = discounted life-years x GDP per capita       (constant 2022 USD, panel)
#   cost per LY      = discounted PEPFAR disbursements / discounted life-years
#                      (pepfar, constant 2022 USD, panel; all PEPFAR spending, not only
#                      maternal health, so this overstates the cost per maternal life-year)
# Every total is linear in (γ1, γ2, γ4), so its standard error comes straight from the
# clustered vcov: se = sqrt(w' V w). X_i = bl_hiv_prev >= 0, so zero means no HIV burden.
#
# Run from the repository root:  Rscript analysis/parallel_trends/life_years_maternal.R
# Downloads two WPP 2024 files into data/raw (or PEPFAR_RAW_DIR) the first time; data/raw
# is gitignored. Writes to tables/cost_effectiveness/{tex, figures, csv}/

source("analysis/econ_setup.R")

# =====================================================
# 0. Settings
# =====================================================
y          <- "mort_maternal"
start_year <- 1990
break_main <- 2008         # main break year
spec_main  <- "Adjusted"   # main specification
age_death  <- 25           # remaining female life expectancy at this exact age (multiple of 5)
r_main     <- 0.03         # discount rate
last_year  <- 2022

raw_dir <- Sys.getenv("PEPFAR_RAW_DIR", "data/raw")
out_dir <- "tables/cost_effectiveness"
tex_dir <- file.path(out_dir, "tex"); fig_dir <- file.path(out_dir, "figures")
csv_dir <- file.path(out_dir, "csv")
for (dd in c(raw_dir, tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Trend-break model (same as trend_break_maternal_2008.R, HIV burden version)
# =====================================================
base_data <- samples[["All PEPFAR"]]   # all PEPFAR recipients + controls; treat = 1 for PEPFAR

add_terms <- function(d, B) {
  d |>
    filter(!is.na(bl_hiv_prev)) |>
    mutate(x         = bl_hiv_prev,
           post_b    = as.integer(year >= B),
           tr_base   = x * (year - start_year),
           tr_change = x * (year - B) * post_b,
           jump      = x * post_b)
}
fe_part <- function(spec) {
  if (spec == "Unadjusted") return("iso3c + year")
  paste0("iso3c + bl_wb_income^year + year[",
         paste(setdiff(bl_cont, "bl_hiv_prev"), collapse = ", "), "]")
}
fit_cache <- list()
fit_model <- function(B, spec) {
  key <- paste(B, spec)
  if (is.null(fit_cache[[key]])) {
    d <- add_terms(base_data, B)
    m <- feols(as.formula(paste0(y, " ~ tr_base + tr_change + jump | ", fe_part(spec))),
               data = d, cluster = ~iso3c)
    fit_cache[[key]] <<- list(m = m, used = d[obs(m), ])
  }
  fit_cache[[key]]
}

# =====================================================
# 2. WPP 2024: births and female life expectancy
# =====================================================
wpp_base <- "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/"
di_file  <- file.path(raw_dir, "WPP2024_Demographic_Indicators_Medium.csv.gz")      # same file build_pepfar_panel_FIXED.R uses
lt_file  <- file.path(raw_dir, "WPP2024_Life_Table_Abridged_Medium_1950-2023.csv.gz")

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

if (!get_file(di_file)) stop("Could not download ", basename(di_file),
                             ". Download it from the WPP 2024 CSV page and save it to ", raw_dir)
di_need <- c("ISO3_code", "Time", "Variant", "Births", "LE15Female", "Q1550Female")
di_miss <- setdiff(di_need, names(read_csv(di_file, n_max = 0, show_col_types = FALSE)))
if (length(di_miss)) stop("Columns missing from the WPP file: ", paste(di_miss, collapse = ", "))
di <- read_csv(di_file, col_select = all_of(di_need), show_col_types = FALSE) |>
  filter(!is.na(ISO3_code), Variant == "Medium", Time >= 2000, Time <= last_year) |>
  transmute(iso3c = ISO3_code, year = as.integer(Time),
            births = Births * 1000,                       # WPP reports thousands
            le15f  = LE15Female,
            q1550f = ifelse(Q1550Female > 1, Q1550Female / 1000, Q1550Female)) |>  # per 1,000 -> share
  distinct(iso3c, year, .keep_all = TRUE)

# Remaining female life expectancy at exact age `age`
lt_need <- c("ISO3_code", "Time", "Sex", "AgeGrpStart", "ex")
lt_ok <- get_file(lt_file) &&
  tryCatch(all(lt_need %in% names(read_csv(lt_file, n_max = 0, show_col_types = FALSE))),
           error = function(e) FALSE)
if (lt_ok) {
  lt <- read_csv(lt_file, col_select = all_of(lt_need), show_col_types = FALSE) |>
    filter(!is.na(ISO3_code), Sex == "Female", Time >= 2000, Time <= last_year) |>
    transmute(iso3c = ISO3_code, year = as.integer(Time), age = AgeGrpStart, ex) |>
    distinct(iso3c, year, age, .keep_all = TRUE)
  ex_source <- "WPP 2024 abridged life table (female)"
  make_ex <- function(age_x) {
    if (!age_x %in% lt$age) stop("age_death must be one of: ", paste(sort(unique(lt$age)), collapse = ", "))
    lt |> filter(age == age_x) |> select(iso3c, year, ex)
  }
} else {
  # Fallback if the life-table file is unavailable: from female e15 and the probability of
  # dying 15-50, assuming a constant death rate between 15 and 50.
  message("Life-table file unavailable; approximating remaining life expectancy from LE15Female.")
  ex_source <- "approximated from WPP 2024 LE15Female and Q1550Female"
  make_ex <- function(age_x) {
    di |> mutate(p  = (1 - q1550f)^((age_x - 15) / 35),
                 ex = (le15f - (age_x - 15) * (1 + p) / 2) / p) |>
      select(iso3c, year, ex)
  }
}
births <- select(di, iso3c, year, births)

# =====================================================
# 3. Counterfactual and totals
# =====================================================

# Country-year rows with every input; k = years since the break.
# scope "sample": countries in the model's estimation sample (the adjusted model drops
# countries with missing baseline covariates, e.g. Tanzania, Uganda, Zambia, Zimbabwe).
# scope "all": also applies the estimated coefficients to those dropped countries.
cy_rows <- function(B, spec, age_x, r, ex_adj, scope = "sample") {
  src <- if (scope == "all") filter(add_terms(base_data, B), !is.na(.data[[y]])) else fit_model(B, spec)$used
  src |>
    filter(year >= B) |>
    select(iso3c, country, year, treat, x, all_of(y), gdp_pc, pepfar) |>
    left_join(births, by = c("iso3c", "year")) |>
    left_join(make_ex(age_x), by = c("iso3c", "year")) |>
    mutate(k      = year - B,
           ex_use = ex * ex_adj,
           ann    = if (r > 0) (1 - (1 + r)^(-ex_use)) / r else ex_use,  # discounted years per death averted
           disc   = (1 + r)^(-k),                                       # back to the break year
           pepfar = coalesce(pepfar, 0))
}

# Weights on (γ1, γ2, γ4) that turn the coefficients into "deaths averted" in each row
row_w <- function(d, cf, jump) {
  base <- d$x * d$births / 1e5
  cbind(tr_base   = if (cf == "flat") -base * d$k else 0,
        tr_change = -base * d$k,
        jump      = if (jump) -base else 0)
}

# Total of f x (deaths averted) over the rows, with its clustered SE
lin_total <- function(W, f, b, V) {
  w <- colSums(W * f)
  nm <- names(w)
  c(est = sum(w * b[nm]), se = sqrt(drop(t(w) %*% V[nm, nm] %*% w)))
}

run_scenario <- function(label, B = break_main, spec = spec_main, age_x = age_death,
                         r = r_main, ex_adj = 1, cf = "trend", jump = TRUE, scope = "sample") {
  fm <- fit_model(B, spec); b <- coef(fm$m); V <- vcov(fm$m)
  all_rows <- cy_rows(B, spec, age_x, r, ex_adj, scope)
  ok <- complete.cases(all_rows[, c("births", "ex", "gdp_pc")])
  d  <- all_rows[ok, ]
  bind_rows(lapply(c("All countries", "PEPFAR recipients"), function(g) {
    dd <- if (g == "PEPFAR recipients") filter(d, treat == 1) else d
    W  <- row_w(dd, cf, jump)
    deaths <- lin_total(W, 1, b, V)
    ly_u   <- lin_total(W, dd$ex_use, b, V)
    ly_d   <- lin_total(W, dd$ann * dd$disc, b, V)
    val    <- lin_total(W, dd$ann * dd$disc * dd$gdp_pc, b, V)
    pep    <- sum(dd$pepfar * dd$disc)
    tibble(scenario = label, group = g, break_year = B, spec = spec, cf = cf,
           jump = jump, scope = scope, n_countries = n_distinct(dd$iso3c), age_death = age_x, r = r, ex_adj = ex_adj,
           country_years = nrow(dd), dropped_missing = sum(!ok),
           obs_deaths     = sum(dd[[y]] * dd$births / 1e5),
           deaths_averted = deaths[["est"]], deaths_se = deaths[["se"]],
           ly_undisc      = ly_u[["est"]],   ly_undisc_se = ly_u[["se"]],
           ly_disc        = ly_d[["est"]],   ly_disc_se   = ly_d[["se"]],
           value_usd      = val[["est"]],    value_se     = val[["se"]],
           pepfar_disc_usd = pep) |>
      mutate(deaths_lo = deaths_averted - 1.96 * deaths_se, deaths_hi = deaths_averted + 1.96 * deaths_se,
             ly_lo     = ly_disc - 1.96 * ly_disc_se,       ly_hi     = ly_disc + 1.96 * ly_disc_se,
             cost_per_ly    = ifelse(ly_disc > 0, pepfar_disc_usd / ly_disc, NA_real_),
             cost_per_ly_lo = ifelse(ly_hi > 0, pepfar_disc_usd / ly_hi, NA_real_),   # more LY -> cheaper
             cost_per_ly_hi = ifelse(ly_lo > 0, pepfar_disc_usd / ly_lo, NA_real_),
             benefit_cost   = value_usd / pepfar_disc_usd)
  }))
}

results <- bind_rows(
  run_scenario("Main"),
  run_scenario("Slope change only (no 2008 jump)", jump = FALSE),
  run_scenario("Flat counterfactual (pre-trend stops at break)", cf = "flat"),
  run_scenario("Undiscounted", r = 0),
  run_scenario("Age at death 20", age_x = 20),
  run_scenario("Age at death 30", age_x = 30),
  run_scenario("Remaining LE x 0.8 (HIV+ women)", ex_adj = 0.8),
  run_scenario("Applied to all recipients (incl. dropped)", scope = "all"),
  run_scenario("Unadjusted model", spec = "Unadjusted"),
  run_scenario("Break in 2004", B = 2004),
  run_scenario("Break in 2012", B = 2012)
)
write_csv(results, file.path(csv_dir, "maternal_life_years_scenarios.csv"))

# =====================================================
# 4. Main scenario by country and by year
# =====================================================
fm <- fit_model(break_main, spec_main); b <- coef(fm$m); V <- vcov(fm$m)
main_rows <- cy_rows(break_main, spec_main, age_death, r_main, 1) |>
  filter(complete.cases(births, ex, gdp_pc))
Wm <- row_w(main_rows, "trend", TRUE)
main_rows <- main_rows |>
  mutate(deaths_obs     = .data[[y]] * births / 1e5,
         deaths_averted = drop(Wm %*% b[colnames(Wm)]),
         deaths_cf      = deaths_obs + deaths_averted,
         ly_disc        = deaths_averted * ann * disc,
         value_usd      = ly_disc * gdp_pc)

by_country <- main_rows |>
  group_by(iso3c, country, treat, hiv_prev_2003 = x) |>
  summarise(deaths_obs = sum(deaths_obs), deaths_averted = sum(deaths_averted),
            ly_disc = sum(ly_disc), value_usd = sum(value_usd),
            pepfar_disc_usd = sum(pepfar * disc), .groups = "drop") |>
  arrange(desc(deaths_averted))
write_csv(by_country, file.path(csv_dir, "maternal_life_years_by_country.csv"))

by_year <- bind_rows(lapply(sort(unique(main_rows$year)), function(t) {
  dd <- filter(main_rows, year == t, treat == 1)
  est <- lin_total(row_w(dd, "trend", TRUE), 1, b, V)
  tibble(year = t, deaths_obs = sum(dd$deaths_obs),
         deaths_averted = est[["est"]], deaths_averted_se = est[["se"]])
})) |>
  mutate(deaths_cf = deaths_obs + deaths_averted,
         cf_lo = deaths_obs + deaths_averted - 1.96 * deaths_averted_se,
         cf_hi = deaths_obs + deaths_averted + 1.96 * deaths_averted_se)
write_csv(by_year, file.path(csv_dir, "maternal_life_years_by_year_pepfar.csv"))

# =====================================================
# 5. Figures
# =====================================================
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}

p_cf <- ggplot(by_year, aes(year)) +
  geom_ribbon(aes(ymin = cf_lo, ymax = cf_hi), fill = pal[2], alpha = 0.2) +
  geom_line(aes(y = deaths_cf, color = "Counterfactual (no break)"), linewidth = 0.9) +
  geom_line(aes(y = deaths_obs, color = "Observed"), linewidth = 0.9) +
  scale_color_manual(values = c("Observed" = pal[1], "Counterfactual (no break)" = pal[2])) +
  scale_y_continuous(labels = function(v) formatC(v, format = "d", big.mark = ",")) +
  labs(title = "Maternal deaths in PEPFAR recipient countries: observed vs. no-break counterfactual",
       subtitle = paste0("Trend break in ", break_main, ", HIV prevalence 2003 x trend, ",
                         tolower(spec_main), " model. Band: 95% CI."),
       x = NULL, y = "Maternal deaths per year") +
  theme_es
save_both(p_cf, "maternal_counterfactual_deaths", 8, 4.5)

top <- by_country |> filter(treat == 1) |> slice_max(deaths_averted, n = 15)
p_top <- ggplot(top, aes(deaths_averted, reorder(country, deaths_averted))) +
  geom_col(fill = pal[1]) +
  scale_x_continuous(labels = function(v) formatC(v, format = "d", big.mark = ",")) +
  labs(title = paste0("Maternal deaths averted, ", break_main, "-", last_year, ": top 15 PEPFAR countries"),
       subtitle = "Point estimates, main scenario", x = "Deaths averted", y = NULL) +
  theme_es + theme(panel.grid.major.x = element_line(color = grid_col, linewidth = 0.3))
save_both(p_top, "maternal_deaths_averted_top15", 7, 5)

# =====================================================
# 6. LaTeX table (PEPFAR recipients)
# =====================================================
fmt0 <- function(v) formatC(round(v), format = "d", big.mark = ",")
fmt_ci <- function(lo, hi) sprintf("[%s, %s]", fmt0(lo), fmt0(hi))
fmt_bn <- function(v) sprintf("%.2f", v / 1e9)
pr <- filter(results, group == "PEPFAR recipients")
rows <- lapply(seq_len(nrow(pr)), function(i) {
  q <- pr[i, ]
  c(q$scenario, fmt0(q$deaths_averted), fmt_ci(q$deaths_lo, q$deaths_hi),
    fmt0(q$ly_disc), fmt_ci(q$ly_lo, q$ly_hi), fmt_bn(q$value_usd),
    ifelse(is.na(q$cost_per_ly), "--", fmt0(q$cost_per_ly)), sprintf("%.2f", q$benefit_cost))
})
write_tex_table(
  file    = file.path(tex_dir, "maternal_life_years.tex"),
  caption = "Maternal deaths and life-years associated with the post-break trend change, PEPFAR recipients",
  label   = "tab:maternal_life_years",
  head    = list(c("", "Deaths", "95\\% CI", "Life-years", "95\\% CI", "Value", "PEPFAR \\$",
                   "Benefit/"),
                 c("Scenario", "averted", "", "(discounted)", "", "(\\$ bn)", "per life-year", "cost")),
  rows    = rows,
  notes   = paste(
    "Counterfactual: each country stays on its pre-break path ($\\gamma_2 = \\gamma_4 = 0$) in the",
    "trend-break model of maternal mortality on HIV prevalence 2003; the flat scenario instead stops",
    "the pre-break trend at the break. Deaths averted = gap in maternal deaths per 100,000 live",
    "births $\\times$ births (UN WPP 2024) / 100,000, summed over PEPFAR recipients from the break",
    "year to", last_year, ". Life-years = deaths averted $\\times$ remaining female life expectancy",
    paste0("at age ", age_death, " (", ex_source, "), discounted at ", 100 * r_main, "\\% unless noted."),
    "Value = discounted life-years $\\times$ GDP per capita; PEPFAR \\$ = discounted PEPFAR",
    "disbursements over the same country-years; both in constant 2022 USD. PEPFAR spending covers",
    "all program areas, so cost per life-year is an upper bound for maternal outcomes alone.",
    "CIs from country-clustered standard errors."),
  align   = "lccccccc")

# =====================================================
# 7. Console summary
# =====================================================
cat("\nLife expectancy source:", ex_source, "\n")
cat("\nPEPFAR recipients, all scenarios:\n")
print(pr |> transmute(scenario, n_countries, deaths_averted = round(deaths_averted), deaths_se = round(deaths_se),
                      ly_disc = round(ly_disc), value_bn = round(value_usd / 1e9, 2),
                      pepfar_bn = round(pepfar_disc_usd / 1e9, 2), cost_per_ly = round(cost_per_ly),
                      benefit_cost = round(benefit_cost, 2)), n = Inf, width = Inf)
cat("\nMain scenario, observed maternal deaths in PEPFAR recipients", break_main, "-", last_year, ":",
    fmt0(pr$obs_deaths[pr$scenario == "Main"]), "\n")
cat("\nDone. Output in", out_dir, "\n")
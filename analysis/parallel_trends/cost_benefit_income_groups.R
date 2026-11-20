# Step 5: cost-effectiveness and cost-benefit by World Bank income group, plus a VSL check
# Author: Rishad Vaghaiwalla
#
# Takes each outcome's main counterfactual (life_years_mortality.R) and, instead of valuing a
# life-year at the country's own GDP per capita (Table 20), values it at World Bank income
# thresholds, separately for each income group, as in the assignment notes:
#   upper-middle income  -> the top of the upper-middle band
#   lower-middle income  -> the top of the lower-middle band
#   low income           -> the low-income threshold
# High-income PEPFAR recipients (a handful) are valued at the upper-middle ceiling.
#
# Panel A, by income group (each country classified in its break year by default):
#   deaths averted (95% CI), discounted life-years, matched PEPFAR dollars,
#   cost-effectiveness = PEPFAR $ per death averted and per life-year,
#   benefit/cost = life-years x threshold / PEPFAR $
# Panel B, all recipients, alternative values of a life-year:
#   own GDP per capita (Table 20), own-group threshold, each threshold applied to everyone,
#   and VSL-based values (Robinson, Hammitt and O'Keeffe 2019, J. Benefit-Cost Analysis):
#     (a) VSL/GNI = 160 x (GNI_i / GNI_US)^0.5 (income elasticity 1.5), floor 20 x GNI_i
#     (b) VSL = 160 x GNI_i       (c) VSL = 100 x GNI_i         [GNI per capita at PPP]
#   VSLY = VSL / remaining life expectancy at age = half of life expectancy at birth
#   (undiscounted, their recommendation), then x discounted life-years like every other row.
#
# Inputs: life_years_setup.R (and Step 4's best_break_year.csv); for the threshold deflation
# and the VSL rows, data/intermediate/us_gdp_deflator.csv and wdi_valuation_2000_2022.csv
# from data_pull/pull_wdi_valuation.R. Without them the VSL rows are skipped and thresholds
# stay in their own year's dollars (a message says so).
#
# Run from the repository root:  Rscript analysis/parallel_trends/cost_benefit_income_groups.R
# Writes to tables/cost_effectiveness/by_outcome/{tex, csv}/

source("analysis/parallel_trends/life_years_setup.R")

outcomes_run <- c("mort_maternal", "mort_child", "mort_allcause")
tex_dir <- file.path(ce_dir, "tex"); csv_dir <- file.path(ce_dir, "csv")
for (dd in c(tex_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 0. Thresholds and settings
# =====================================================
# World Bank FY2027 classification (July 2026), Atlas GNI per capita of 2025, current US$:
# low <= 1,175; lower-middle 1,176-4,635; upper-middle 4,636-14,375; high > 14,375.
thr_raw  <- c(L = 1175, LM = 4635, UM = 14375)
thr_year <- 2025              # dollars the thresholds are in
thr_label <- "World Bank FY2027 thresholds"
class_year <- "break"         # income group as of: "break" (break year), "2022", or "each" (year by year)

defl_file <- "data/intermediate/us_gdp_deflator.csv"
wdi_file  <- "data/intermediate/wdi_valuation_2000_2022.csv"
defl <- if (file.exists(defl_file)) read_csv(defl_file, show_col_types = FALSE) else NULL
wdi  <- if (file.exists(wdi_file)) read_csv(wdi_file, show_col_types = FALSE) else NULL

to_2022 <- function(v, yr) v * defl$deflator[defl$year == 2022] / defl$deflator[match(yr, defl$year)]
if (!is.null(defl)) {
  use_yr <- if (thr_year %in% defl$year) thr_year else max(defl$year)
  if (use_yr != thr_year) message("US deflator stops in ", use_yr, "; deflating the thresholds from ", use_yr)
  thr <- to_2022(thr_raw, use_yr)
  thr_note <- paste0(thr_label, " (", thr_year, " US\\$) converted to 2022 US\\$ with the US GDP deflator")
} else {
  message("Note: ", defl_file, " not found (run data_pull/pull_wdi_valuation.R); thresholds left in ",
          thr_year, " dollars.")
  thr <- thr_raw
  thr_note <- paste0(thr_label, " in ", thr_year, " US\\$ (not deflated)")
}
grp_lab <- c(L = "Low income", LM = "Lower-middle income", UM = "Upper-middle income (and high)")

# =====================================================
# 1. VSL per life-year (Robinson, Hammitt and O'Keeffe 2019), per country-year
# =====================================================
ex_total <- lt |> filter(sex == "Total") |> select(iso3c, year, age, ex)
vsly_cols <- function(d) {
  if (is.null(wdi) || is.null(defl)) return(d)
  us <- wdi |> filter(iso3c == "USA") |> select(year, gni_us = gni_pc_ppp)
  d |>
    left_join(select(wdi, iso3c, year, gni_pc_ppp), by = c("iso3c", "year")) |>
    left_join(us, by = "year") |>
    mutate(age_star = 5 * round(le0 / 2 / 5)) |>                         # half of life expectancy at birth
    left_join(rename(ex_total, ex_star = ex), by = c("iso3c", "year", "age_star" = "age")) |>
    mutate(ratio_a = pmax(20, 160 * (gni_pc_ppp / gni_us)^0.5),
           vsl_a = ratio_a * gni_pc_ppp, vsl_b = 160 * gni_pc_ppp, vsl_c = 100 * gni_pc_ppp,
           across(c(vsl_a, vsl_b, vsl_c), ~ to_2022(.x, year)),         # to 2022 dollars
           vsly_a = vsl_a / ex_star, vsly_b = vsl_b / ex_star, vsly_c = vsl_c / ex_star)
}

# =====================================================
# 2. Totals by income group and by valuation, one outcome
# =====================================================
run_outcome <- function(y) {
  mr <- main_rows(y)
  d  <- mutate(mr$rows, .i = row_number())
  grp_of <- function(g) recode(as.character(g), H = "UM")
  if (class_year == "each") {
    d$grp <- grp_of(d$wb_income)
  } else {
    cy <- if (class_year == "break") mr$B else as.integer(class_year)
    cls <- panel |> filter(year == cy) |> select(iso3c, grp0 = wb_income)
    d <- d |> left_join(cls, by = "iso3c") |> mutate(grp = grp_of(grp0)) |> select(-grp0)
  }
  d <- vsly_cols(d) |> filter(treat == 1)
  if (any(is.na(d$grp))) message(y, ": ", n_distinct(d$iso3c[is.na(d$grp)]),
                                 " recipient(s) with no income group are left out of Panel A")
  W <- mr$W[d$.i, , drop = FALSE]
  tot <- function(rows, value, label, panel_name) {
    rows <- rows & !is.na(value)
    totals(d[rows, ], W[rows, , drop = FALSE], mr$b, mr$V, value[rows]) |>
      mutate(outcome = y, break_year = mr$B, panel = panel_name, row = label,
             value_per_ly = sum(d$ann[rows] * d$disc[rows] * value[rows] * d$deaths_averted[rows]) /
                            sum(d$ann[rows] * d$disc[rows] * d$deaths_averted[rows]),   # life-year-weighted mean
             .before = 1)
  }
  thr_row <- unname(thr[d$grp])
  A <- bind_rows(lapply(c("L", "LM", "UM"), function(g)
         tot(d$grp %in% g, thr_row, grp_lab[[g]], "A")),
       tot(!is.na(d$grp), thr_row, "All recipients", "A"))
  B <- bind_rows(
    tot(rep(TRUE, nrow(d)), d$gdp_pc, "Own GDP per capita (Table 20)", "B"),
    tot(!is.na(d$grp), thr_row, "Own income group's threshold", "B"),
    tot(rep(TRUE, nrow(d)), rep(thr[["L"]], nrow(d)), "Everyone at the low-income threshold", "B"),
    tot(rep(TRUE, nrow(d)), rep(thr[["LM"]], nrow(d)), "Everyone at the lower-middle ceiling", "B"),
    tot(rep(TRUE, nrow(d)), rep(thr[["UM"]], nrow(d)), "Everyone at the upper-middle ceiling", "B"),
    if ("vsly_a" %in% names(d)) bind_rows(
      tot(rep(TRUE, nrow(d)), d$vsly_a, "VSL (a): ratio 160, elasticity 1.5, floor 20", "B"),
      tot(rep(TRUE, nrow(d)), d$vsly_b, "VSL (b): 160 x GNI per capita", "B"),
      tot(rep(TRUE, nrow(d)), d$vsly_c, "VSL (c): 100 x GNI per capita", "B")))
  bind_rows(A, B)
}

res <- bind_rows(lapply(outcomes_run, run_outcome))
write_csv(res, file.path(csv_dir, "cost_benefit_income_groups.csv"))

# =====================================================
# 3. LaTeX, one table per outcome
# =====================================================
fmt_usd <- function(v) ifelse(is.na(v), "--", paste0("\\$", fmt0(v)))
for (y in outcomes_run) {
  cfg <- outcome_cfg[[y]]; r <- filter(res, outcome == y)
  rowsA <- lapply(which(r$panel == "A"), function(i) { q <- r[i, ]
    c(q$row, as.character(q$n_countries), fmt0(q$deaths_averted), fmt_ci(q$deaths_lo, q$deaths_hi),
      fmt0(q$ly_disc), fmt_bn(q$pepfar_match_disc_usd), fmt_usd(q$cost_per_death_match),
      fmt_usd(q$cost_per_ly_match), fmt_usd(q$value_per_ly), sprintf("%.2f", q$benefit_cost_match)) })
  rowsB <- lapply(which(r$panel == "B"), function(i) { q <- r[i, ]
    c(q$row, "", "", "", "", "", "", "", fmt_usd(q$value_per_ly), sprintf("%.2f", q$benefit_cost_match)) })
  write_tex_table(
    file    = file.path(tex_dir, paste0(cfg$file, "_cost_benefit.tex")),
    caption = paste0(cfg$label, " mortality: cost-effectiveness and benefit/cost by World Bank income group"),
    label   = paste0("tab:", cfg$file, "_cost_benefit"),
    head    = list(c("", "", "Deaths", "", "Life-years", "PEPFAR \\$", "\\$ per death", "\\$ per",
                     "Value per", "Benefit/"),
                   c("", "Countries", "averted", "95\\% CI", "(disc.)", "(bn)", "averted", "life-year",
                     "life-year", "cost")),
    rows    = c(list("\\multicolumn{10}{l}{\\textit{A. By income group, life-year valued at the group's threshold}} \\\\"),
                rowsA, list("\\midrule",
                "\\multicolumn{10}{l}{\\textit{B. All recipients, alternative values of a life-year}} \\\\"), rowsB),
    notes   = paste(
      paste0("Main scenario of the ", tolower(cfg$label), " life-years table (break in ", r$break_year[1],
             ", HIV prevalence 2003 $\\times$ trend, adjusted, no-break counterfactual, 3\\% discounting),"),
      "PEPFAR recipients. Income group as of", ifelse(class_year == "break", "the break year",
                                                     ifelse(class_year == "each", "each year", class_year)),
      paste0("(World Bank historical classification); high-income recipients are valued at the upper-middle ceiling."),
      paste0("Thresholds: ", thr_note, ". PEPFAR \\$ = discounted ", dollars_note(y), ", constant 2022 US\\$."),
      "Cost-effectiveness = PEPFAR \\$ per death averted and per discounted life-year. Benefit/cost =",
      "discounted life-years $\\times$ value per life-year / PEPFAR \\$. Value per life-year in Panel B is the",
      "life-year-weighted mean. VSL rows follow the standardized sensitivity analysis of Robinson, Hammitt",
      "and O'Keeffe (2019) with GNI per capita at PPP, converted to a value per life-year by dividing by",
      "remaining life expectancy at half of life expectancy at birth."),
    align   = "lccccccccc")
}

# =====================================================
# 4. Console summary
# =====================================================
cat("\nThresholds used (per life-year):", paste(names(thr), round(thr), sep = " = ", collapse = ", "), "\n")
for (y in outcomes_run) {
  cat("\n==========", outcome_cfg[[y]]$label, "==========\n")
  print(filter(res, outcome == y) |>
          transmute(panel, row, n_countries, deaths = round(deaths_averted), ly_m = round(ly_disc / 1e6, 2),
                    pepfar_bn = round(pepfar_match_disc_usd / 1e9, 2), usd_per_death = round(cost_per_death_match),
                    usd_per_ly = round(cost_per_ly_match), value_per_ly = round(value_per_ly),
                    bc = round(benefit_cost_match, 2)), n = Inf, width = Inf)
}
cat("\nDone. Output in", ce_dir, "\n")

# Step 2: Table 20 and Figure 5 for child, all-cause and maternal mortality
# Author: Rishad Vaghaiwalla
#
# life_years_maternal.R, run for every mortality outcome. Same model (HIV prevalence 2003 x
# trend, adjusted, country-clustered SEs), same counterfactual (no break: γ2 = γ4 = 0), same
# 11 scenarios, same discounting (3%) and dollars (constant 2022 USD). Per outcome only the
# rate, its denominator and the life expectancy change (life_years_setup.R, outcome_cfg):
#
#   Maternal   per 100,000 live births x WPP births   female life expectancy at 25   break 2008
#   Child      per 1,000 live births   x WPP births   life expectancy at birth        break: Step 4
#   All-cause  per 1,000 people        x WPP pop.     life expectancy at 35           break: Step 4
#
# Maternal is included as a check: with no shares file its numbers equal Table 20.
# With data/intermediate/pepfar_outcome_shares.csv (Step 3), maternal and child costs use
# matched PEPFAR dollars; all-cause always uses all PEPFAR dollars.
#
# Run order: trend_break_sweep_mortality.R (Step 4) first, then this script.
# Run from the repository root:  Rscript analysis/parallel_trends/life_years_mortality.R
# Writes to tables/cost_effectiveness/by_outcome/{tex, figures, csv}/

source("analysis/parallel_trends/life_years_setup.R")

outcomes_run <- c("mort_child", "mort_allcause", "mort_maternal")
tex_dir <- file.path(ce_dir, "tex"); fig_dir <- file.path(ce_dir, "figures"); csv_dir <- file.path(ce_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. The 11 Table 20 scenarios, for one outcome
# =====================================================
scenarios <- function(y) {
  cfg <- outcome_cfg[[y]]; B <- main_break(y)
  # Two alternative break years from {2004, 2008, 2012}: the nearest below and above the main
  # year (maternal: 2004 and 2012, as in Table 20), never before PEPFAR started in 2004
  cand  <- setdiff(c(2004, 2008, 2012), B)
  below <- cand[cand < B]; above <- cand[cand > B]
  alt <- if (length(below) && length(above)) c(max(below), min(above)) else head(sort(cand), 2)
  bind_rows(
    run_scenario(y, "Main"),
    run_scenario(y, sprintf("Slope change only (no %d jump)", B), jump = FALSE),
    run_scenario(y, "Flat counterfactual (pre-trend stops at break)", cf = "flat"),
    run_scenario(y, "Undiscounted", r = 0),
    run_scenario(y, sprintf("Age at death %d", cfg$ages_alt[1]), age_x = cfg$ages_alt[1]),
    run_scenario(y, sprintf("Age at death %d", cfg$ages_alt[2]), age_x = cfg$ages_alt[2]),
    run_scenario(y, "Remaining LE x 0.8 (HIV+)", ex_adj = 0.8),
    run_scenario(y, "Applied to all recipients (incl. dropped)", scope = "all"),
    run_scenario(y, "Unadjusted model", spec = "Unadjusted"),
    run_scenario(y, sprintf("Break in %d", alt[1]), B = alt[1]),
    run_scenario(y, sprintf("Break in %d", alt[2]), B = alt[2])
  )
}

all_results <- list()
for (y in outcomes_run) {
  cfg <- outcome_cfg[[y]]; B <- main_break(y); f <- cfg$file
  cat("\n==========", cfg$label, "mortality, break in", B, "==========\n")

  results <- scenarios(y)
  write_csv(results, file.path(csv_dir, paste0(f, "_life_years_scenarios.csv")))
  all_results[[y]] <- results

  # =====================================================
  # 2. Main scenario by country and by year (from life_years_maternal.R, section 4)
  # =====================================================
  mr <- main_rows(y)
  by_country <- mr$rows |>
    group_by(iso3c, country, treat, hiv_prev_2003 = x) |>
    summarise(deaths_obs = sum(deaths_obs), deaths_averted = sum(deaths_averted),
              ly_disc = sum(ly_disc), value_usd = sum(value_usd),
              pepfar_disc_usd = sum(pepfar * disc), pepfar_match_disc_usd = sum(pepfar_match * disc),
              .groups = "drop") |>
    arrange(desc(deaths_averted))
  write_csv(by_country, file.path(csv_dir, paste0(f, "_life_years_by_country.csv")))

  by_year <- by_year_cf(mr)
  write_csv(by_year, file.path(csv_dir, paste0(f, "_life_years_by_year_pepfar.csv")))

  # =====================================================
  # 3. Figures (Figure 5 and the top-15 bar chart)
  # =====================================================
  p_cf <- ggplot(by_year, aes(year)) +
    geom_ribbon(aes(ymin = cf_lo, ymax = cf_hi), fill = pal[2], alpha = 0.2) +
    geom_line(aes(y = deaths_cf, color = "Counterfactual (no break)"), linewidth = 0.9) +
    geom_line(aes(y = deaths_obs, color = "Observed"), linewidth = 0.9) +
    scale_color_manual(values = c("Observed" = pal[1], "Counterfactual (no break)" = pal[2])) +
    scale_y_continuous(labels = big_lab) +
    labs(title = paste0(cfg$label, " deaths, PEPFAR recipients: observed vs. no-break counterfactual"),
         subtitle = paste0("Trend break in ", B, ", HIV prevalence 2003 x trend, ",
                           tolower(spec_main), " model. Band: 95% CI."),
         x = NULL, y = paste0(cfg$label, " deaths per year")) +
    theme_es
  save_both(p_cf, fig_dir, paste0(f, "_counterfactual_deaths"), 8, 4.5)

  top <- by_country |> filter(treat == 1) |> slice_max(deaths_averted, n = 15)
  p_top <- ggplot(top, aes(deaths_averted, reorder(country, deaths_averted))) +
    geom_col(fill = pal[1]) +
    scale_x_continuous(labels = big_lab) +
    labs(title = paste0(cfg$label, " deaths averted, ", B, "-", last_year, ": top 15 PEPFAR countries"),
         subtitle = "Point estimates, main scenario", x = "Deaths averted", y = NULL) +
    theme_es + theme(panel.grid.major.x = element_line(color = grid_col, linewidth = 0.3))
  save_both(p_top, fig_dir, paste0(f, "_deaths_averted_top15"), 7, 5)

  # =====================================================
  # 4. LaTeX table (PEPFAR recipients), same layout as Table 20
  # =====================================================
  pr <- filter(results, group == "PEPFAR recipients")
  rows <- lapply(seq_len(nrow(pr)), function(i) {
    q <- pr[i, ]
    c(q$scenario, fmt0(q$deaths_averted), fmt_ci(q$deaths_lo, q$deaths_hi),
      fmt0(q$ly_disc), fmt_ci(q$ly_lo, q$ly_hi), fmt_bn(q$value_usd),
      ifelse(is.na(q$cost_per_ly_match), "--", fmt0(q$cost_per_ly_match)),
      sprintf("%.2f", q$benefit_cost_match))
  })
  write_tex_table(
    file    = file.path(tex_dir, paste0(f, "_life_years.tex")),
    caption = paste0(cfg$label, " deaths and life-years associated with the post-break trend change, PEPFAR recipients"),
    label   = paste0("tab:", f, "_life_years"),
    head    = list(c("", "Deaths", "95\\% CI", "Life-years", "95\\% CI", "Value", "PEPFAR \\$", "Benefit/"),
                   c("Scenario", "averted", "", "(discounted)", "", "(\\$ bn)", "per life-year", "cost")),
    rows    = rows,
    notes   = paste(
      "Counterfactual: each country stays on its pre-break path ($\\gamma_2 = \\gamma_4 = 0$) in the",
      paste0("trend-break model (break in ", B, ") of ", cfg$rate_tex, " on HIV prevalence 2003;"),
      "the flat scenario instead stops the pre-break trend at the break. Deaths averted = gap in the",
      paste0("rate $\\times$ ", ifelse(cfg$denom == "births", "births", "population"), " (UN WPP 2024) /"),
      paste0(fmt0(cfg$scale), ", summed over PEPFAR recipients from the break year to ", last_year, "."),
      paste0("Life-years = deaths averted $\\times$ ", cfg$le_note, " at age ", cfg$age, " (", ex_source, "),"),
      paste0("discounted at ", 100 * r_main, "\\% unless noted. Value = discounted life-years $\\times$ GDP per capita;"),
      paste0("PEPFAR \\$ = discounted ", dollars_note(y), " over the same country-years;"),
      "both in constant 2022 USD. CIs from country-clustered standard errors."),
    align   = "lccccccc")

  cat("\nPEPFAR recipients, all scenarios:\n")
  print(pr |> transmute(scenario, n_countries, deaths_averted = round(deaths_averted),
                        deaths_se = round(deaths_se), ly_disc = round(ly_disc),
                        value_bn = round(value_usd / 1e9, 2), pepfar_bn = round(pepfar_match_disc_usd / 1e9, 2),
                        cost_per_death = round(cost_per_death_match), cost_per_ly = round(cost_per_ly_match),
                        benefit_cost = round(benefit_cost_match, 2)), n = Inf, width = Inf)
}

# =====================================================
# 5. One summary table: main scenario, all three outcomes side by side
# =====================================================
main <- bind_rows(all_results) |> filter(scenario == "Main", group == "PEPFAR recipients")
write_csv(main, file.path(csv_dir, "main_scenario_by_outcome.csv"))
# Cost side as a bracket: matched dollars (the program aimed at the outcome; credits it with
# every death averted, so the cost per death is a LOWER bound) and all PEPFAR dollars (an
# UPPER bound). For all-cause mortality the two are the same.
dash <- function(v, f = fmt0) ifelse(is.na(v), "--", f(v))
col_of <- function(y) {
  q <- filter(main, outcome == y)
  c(as.character(q$break_year), fmt0(q$deaths_averted), fmt_ci(q$deaths_lo, q$deaths_hi),
    fmt0(q$ly_disc),
    fmt_bn(q$pepfar_match_disc_usd), dash(q$cost_per_death_match), dash(q$cost_per_ly_match),
    sprintf("%.2f", q$benefit_cost_match),
    fmt_bn(q$pepfar_disc_usd), dash(ifelse(q$deaths_averted > 0, q$pepfar_disc_usd / q$deaths_averted, NA)),
    dash(q$cost_per_ly), sprintf("%.2f", q$benefit_cost))
}
lab <- c("Break year", "Deaths averted", "95\\% CI", "Life-years (discounted)",
         "\\textit{Program dollars aimed at the outcome}", "\\quad PEPFAR \\$ (bn)", "\\quad \\$ per death averted",
         "\\quad \\$ per life-year", "\\quad Benefit/cost",
         "\\textit{All PEPFAR dollars}", "\\quad PEPFAR \\$ (bn)", "\\quad \\$ per death averted",
         "\\quad \\$ per life-year", "\\quad Benefit/cost")
cols <- sapply(outcomes_run, col_of)
cols <- rbind(cols[1:4, ], "", cols[5:8, ], "", cols[9:12, ])
write_tex_table(
  file    = file.path(tex_dir, "main_scenario_by_outcome.tex"),
  caption = "Deaths averted, cost-effectiveness and benefit/cost by mortality outcome, PEPFAR recipients",
  label   = "tab:life_years_by_outcome",
  head    = list(c("", sapply(outcomes_run, function(y) outcome_cfg[[y]]$label))),
  rows    = lapply(seq_along(lab), function(i) c(lab[i], cols[i, ])),
  notes   = paste(
    "Main scenario of each outcome's life-years table (HIV prevalence 2003 $\\times$ trend, adjusted,",
    "no-break counterfactual, 3\\% discounting, constant 2022 USD). Child and all-cause break years are",
    "the best-fitting years from 2004 on (break-year sweep); maternal uses 2008.",
    "Program dollars: PEPFAR disbursements times the share of PEPFAR program spending on PMTCT (maternal)",
    "or on PMTCT, pediatric care and treatment and orphans and vulnerable children (under-5), from PEPFAR",
    "operational plans and expenditure data (data/intermediate/pepfar\\_outcome\\_shares.csv); all-cause uses",
    "all PEPFAR dollars. Crediting all deaths averted to the targeted program makes its cost per death a lower",
    "bound; dividing by all PEPFAR dollars gives an upper bound.",
    if (is.null(shares)) "NOTE: shares file missing, so both rows use all PEPFAR dollars." else ""),
  align = "lccc", resize = TRUE)

cat("\nDone. Output in", ce_dir, "\n")

# Step 2 (notes: "counterfactual for mortality in the regressions so far, plots if applicable")
# Author: Rishad Vaghaiwalla
#
# For every trend-break regression already in trend_break_mortality.R (3 outcomes x
# 3 treatment versions x 2 specifications, break in 2004), plus each outcome's main break
# year, this plots observed deaths in the treated countries against the no-break
# counterfactual (γ2 = γ4 = 0), with a 95% band. Same machinery as life_years_mortality.R.
#
#   All PEPFAR  treated = PEPFAR recipients, X_i = 1
#   COP         treated = COP countries,     X_i = 1
#   HIV burden  treated = PEPFAR recipients, X_i = HIV prevalence 2003
#
# Run from the repository root:  Rscript analysis/parallel_trends/counterfactual_plots_mortality.R
# Writes to tables/cost_effectiveness/by_outcome/{figures, csv}/

source("analysis/parallel_trends/life_years_setup.R")

fig_dir <- file.path(ce_dir, "figures"); csv_dir <- file.path(ce_dir, "csv")
for (dd in c(fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

outcomes_run <- c("mort_maternal", "mort_child", "mort_allcause")
specs <- c("Unadjusted", "Adjusted")

# break_set: "2004" = the regressions in trend_break_mortality.R; "main" = each outcome's
# main break year (maternal 2008, child and all-cause from the Step 4 sweep)
cf_all <- bind_rows(lapply(c("2004", "main"), function(bs)
  bind_rows(lapply(outcomes_run, function(y)
    bind_rows(lapply(names(versions), function(v)
      bind_rows(lapply(specs, function(s) {
        B <- if (bs == "2004") 2004 else main_break(y)
        by_year_cf(main_rows(y, B, s, v)) |>
          mutate(outcome = outcome_cfg[[y]]$label, version = v, spec = s, break_year = B, break_set = bs)
      })))))))) |>
  mutate(outcome = factor(outcome, levels = sapply(outcomes_run, function(y) outcome_cfg[[y]]$label)),
         version = factor(version, levels = names(versions)))
write_csv(cf_all, file.path(csv_dir, "counterfactual_by_year_all_regressions.csv"))

for (bs in c("2004", "main")) for (s in specs) {
  d <- filter(cf_all, break_set == bs, spec == s)
  sub <- if (bs == "2004") "Break in 2004 (the regressions in trend_break_mortality.R)" else
    "Main break year per outcome (maternal 2008; child and all-cause from the break-year sweep)"
  p <- ggplot(d, aes(year)) +
    geom_vline(aes(xintercept = break_year), color = ink2, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(aes(ymin = cf_lo, ymax = cf_hi), fill = pal[2], alpha = 0.2) +
    geom_line(aes(y = deaths_cf, color = "Counterfactual (no break)"), linewidth = 0.7) +
    geom_line(aes(y = deaths_obs, color = "Observed"), linewidth = 0.7) +
    scale_color_manual(values = c("Observed" = pal[1], "Counterfactual (no break)" = pal[2])) +
    scale_y_continuous(labels = function(v) paste0(formatC(v / 1e6, format = "f", digits = 1), "m")) +
    facet_grid(outcome ~ version, scales = "free_y") +
    labs(title = paste0("Observed vs. no-break counterfactual deaths in treated countries, ", tolower(s), " model"),
         subtitle = paste0(sub, ". Band: 95% CI."), x = NULL, y = "Deaths per year (millions)") +
    theme_es
  save_both(p, fig_dir, paste0("counterfactual_grid_break", bs, "_", tolower(s)), 10, 8)
}

cat("\nDeaths averted over the whole post-break period, treated countries (point estimates):\n")
print(cf_all |> group_by(break_set, outcome, version, spec, break_year) |>
        summarise(deaths_averted = round(sum(deaths_averted)), .groups = "drop") |>
        tidyr::pivot_wider(names_from = spec, values_from = deaths_averted), n = Inf, width = Inf)
cat("\nDone. Output in", ce_dir, "\n")

# Runs the CREL cost-effectiveness assignment end to end, in order
# Author: Rishad Vaghaiwalla
#
# From the repository root:  Rscript analysis/run_cost_effectiveness.R
# or in RStudio (working directory = the repo):  source("analysis/run_cost_effectiveness.R")
#
# Optional steps are skipped with a message if their inputs are missing; the rest still run:
#   - pull_wdi_valuation.R needs the internet (World Bank API); without it Step 5 skips the
#     VSL rows and leaves the income thresholds undeflated.
#   - build_pepfar_outcome_shares.R downloads PEPFAR spending data from data.pepfar.gov into
#     data/raw/pepfar_spending/; if it can't, maternal and child costs fall back to all PEPFAR dollars.
steps <- list(
  list(file = "data_pull/pull_wdi_valuation.R",                          what = "Step 5 inputs: GNI per capita (PPP), US deflator", optional = TRUE),
  list(file = "data_pull/build_pepfar_outcome_shares.R",                 what = "Step 3: PEPFAR dollars matched to maternal and child", optional = TRUE),
  list(file = "analysis/parallel_trends/trend_break_sweep_mortality.R",  what = "Step 4: break-year sweep", optional = FALSE),
  list(file = "analysis/parallel_trends/life_years_mortality.R",         what = "Step 2: Table 20 and Figure 5 for each outcome", optional = FALSE),
  list(file = "analysis/parallel_trends/counterfactual_plots_mortality.R", what = "Step 2: counterfactual plots, all regressions", optional = FALSE),
  list(file = "analysis/parallel_trends/cost_benefit_income_groups.R",   what = "Step 5: CE and benefit/cost by income group", optional = FALSE),
  list(file = "analysis/counterfactual_epi/cost_to_end_hiv.R",           what = "Step 6: cost of ending HIV", optional = FALSE)
)

for (s in steps) {
  message("\n>>> ", s$what, "  (", s$file, ")")
  ok <- tryCatch({ source(s$file, local = new.env()); TRUE },
                 error = function(e) {
                   if (!s$optional) stop(s$file, " failed: ", conditionMessage(e), call. = FALSE)
                   message("    skipped: ", conditionMessage(e)); FALSE
                 })
}
message("\nAll steps done.")

# LaTeX table of the linear-time-trend checks for the Crown economic outcomes
# Author: Anthony Babu
#
# Reads the numbers that event_study_economic.R already wrote (section 3, linear time
# trends) and lays them out as one table per treatment group. Columns:
#   (1) DiD          Y = a_i + g_t + beta*(PEPFAR x Post)                  no trend
#   (2) Trend        theta in Y = a_i + g_t + theta*(PEPFAR x t) + delta*(PEPFAR x Post)
#                    = how much faster PEPFAR countries drift each year, over 1990-2022
#   (3) DiD + trend  delta from the same model: the post-2004 jump on top of that drift
#   (4) Country trends  delta when every country gets its own linear trend
# All unadjusted, SEs clustered by country.
#
# Run from the repository root (after event_study_economic.R):
#   Rscript analysis/parallel_trends/linear_trends_table.R
# Writes to tables/crown_event_studies/robustness/tex/

source("analysis/econ_setup.R")

in_dir  <- "tables/crown_event_studies/robustness"
out_dir <- file.path(in_dir, "tex")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trends <- read_csv(file.path(in_dir, "linear_trends.csv"), show_col_types = FALSE)
summ   <- read_csv(file.path(in_dir, "es_summary.csv"), show_col_types = FALSE) |>
  filter(spec == "Unadjusted") |>
  select(outcome, group, did_twfe, did_twfe_se, pre_wald_p)
tab <- left_join(trends, summ, by = c("outcome", "group"))

# coefficient with stars from est / se (normal p-value), and its SE line
coef_cell <- function(est, se) tex_coef(est, 2 * pnorm(-abs(est / se)))

for (g in groups) {
  rows <- list()
  for (y in names(outcomes)) {
    r <- filter(tab, group == g, outcome == y)
    rows[[length(rows) + 1]] <- c(outcomes_tex[[y]],
      coef_cell(r$did_twfe, r$did_twfe_se), coef_cell(r$slide_theta, r$slide_theta_se),
      coef_cell(r$slide_delta, r$slide_delta_se), coef_cell(r$ctrytrend_delta, r$ctrytrend_delta_se),
      if (r$pre_wald_p < 0.001) "$<$0.001" else sprintf("%.3f", r$pre_wald_p))
    rows[[length(rows) + 1]] <- c("", tex_se(r$did_twfe_se), tex_se(r$slide_theta_se),
      tex_se(r$slide_delta_se), tex_se(r$ctrytrend_delta_se), "")
  }
  write_tex_table(
    file    = file.path(out_dir, paste0("linear_trends_", group_file[[g]], ".tex")),
    caption = paste0("DiD with linear time trends: ", g, " vs.\\ controls"),
    label   = paste0("tab:crown_linear_trends_", group_file[[g]]),
    head    = list(c("", "(1)", "(2)", "(3)", "(4)", ""),
                   c("", "DiD", "PEPFAR trend", "DiD + trend", "Country trends", "Pre-test $p$")),
    rows    = rows,
    notes   = paste(
      "Column (1): two-way fixed effects DiD, $Y_{it} = \\alpha_i + \\gamma_t + \\beta(\\mathrm{PEPFAR}_i \\times \\mathrm{Post}_t)$.",
      "Columns (2)--(3): the same model plus a PEPFAR-specific linear trend,",
      "$Y_{it} = \\alpha_i + \\gamma_t + \\theta(\\mathrm{PEPFAR}_i \\times t) + \\delta(\\mathrm{PEPFAR}_i \\times \\mathrm{Post}_t)$;",
      "(2) reports $\\theta$, the extra change per year in PEPFAR countries over 1990--2022, and (3) reports $\\delta$,",
      "the post-2004 shift on top of that trend. Column (4): $\\delta$ when every country has its own linear trend.",
      "Pre-test $p$: joint F-test that the 13 pre-2004 event-study coefficients are zero (unadjusted).",
      "No covariates. Sample 1990--2022; controls are the 64 comparison countries valid through 2022.",
      "Standard errors clustered by country in parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lccccc")
}
cat("Linear-trend tables written to", out_dir, "\n")

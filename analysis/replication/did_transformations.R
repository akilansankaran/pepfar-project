# Does the Crown DiD survive transforming the outcome? Raw, winsorized, log(x+1), IHS
# Author: Anthony Babu
#
# Same regression as Panel A of did_econ.R (Crown's pooled OLS with all 11 baseline
# covariates), re-run with the outcome transformed four ways:
#   Raw         the outcome as in did_econ.R (levels for GDP per capita and GDP)
#   Winsorized  values below the 1st / above the 99th percentile set to those
#               percentiles (pooled over the whole panel), to blunt outliers such as
#               Equatorial Guinea's +140% growth in 1997
#   log(x+1)    log(1 + y) for levels. Growth rates are in %, and log(1 + y) would drop
#               every year with growth at or below -1% (about a fifth of the panel), so
#               for them we use 100 * log(1 + y/100), the growth rate in log points
#   IHS         inverse hyperbolic sine, asinh(y) = log(y + sqrt(y^2 + 1)); like a log
#               for large y but defined at zero and for negative values
# With log or IHS, the DiD coefficient is roughly a proportional change (0.10 ~ 10%)
# for large values of y; for growth rates near zero, IHS mainly pulls in the extremes.
#
# For each transformation we also re-run the event study (unadjusted and adjusted)
# and test parallel trends: a joint F-test that all 13 pre-2004 coefficients are zero.
#
# Run from the repository root:  Rscript analysis/replication/did_transformations.R
# Writes to tables/crown_did/transformations/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_did/transformations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Outcomes in raw form (GDP per capita and GDP in levels, not logs)
raw_outcomes <- c(
  gdp_pc_growth        = "GDP pc growth",
  gdp_growth           = "GDP growth",
  emp_female           = "Female emp.",
  emp_male             = "Male emp.",
  prim_disengage_girls = "Girls' diseng.",
  prim_disengage_boys  = "Boys' diseng.",
  gdp_pc               = "GDP per capita",
  gdp_bn               = "GDP"
)
growth_rates <- c("gdp_pc_growth", "gdp_growth")

winsor <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

# Build every transformed column once, on the full panel (so the winsorizing
# cut-offs are the same for both treatment groups)
transforms <- c(raw = "Raw", win = "Winsorized", log = "log(x+1)", ihs = "IHS")
panel_t <- panel
for (y in names(raw_outcomes)) {
  panel_t[[paste0(y, "__raw")]] <- panel[[y]]
  panel_t[[paste0(y, "__win")]] <- winsor(panel[[y]])
  panel_t[[paste0(y, "__log")]] <- if (y %in% growth_rates) 100 * log1p(panel[[y]] / 100) else log1p(panel[[y]])
  panel_t[[paste0(y, "__ihs")]] <- asinh(panel[[y]])
}

# Crown's pooled OLS with the 11 baseline covariates (Panel A of did_econ.R)
run_did_crown <- function(d, yvar) {
  f <- as.formula(paste(yvar, "~ post + treat + post:treat +", paste(bl_all, collapse = " + ")))
  feols(f, data = d, cluster = ~iso3c)
}

res <- list()
for (g in groups) {
  d <- make_sample(g, data = panel_t)
  for (y in names(raw_outcomes)) for (tr in names(transforms)) {
    yvar <- paste0(y, "__", tr)
    m <- run_did_crown(d, yvar)
    ct <- coeftable(m)["post:treat", ]
    pre_u <- pre_test(run_es(d, yvar, adjusted = FALSE))
    pre_a <- pre_test(run_es(d, yvar, adjusted = TRUE))
    res[[length(res) + 1]] <- tibble(
      outcome = y, group = g, transform = tr,
      est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
      n_obs = nobs(m), n_countries = n_countries(m, d),
      pre_F_unadj = pre_u[["stat"]], pre_p_unadj = pre_u[["p"]],
      pre_F_adj   = pre_a[["stat"]], pre_p_adj   = pre_a[["p"]])
  }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "did_transformations.csv"))

# ---- LaTeX: one table per group, rows = outcomes, columns = transformations ----
cell <- function(g, y, tr, what) {
  r <- filter(res, group == g, outcome == y, transform == tr)
  if (nrow(r) == 0) return(if (what == "coef") "---" else "")
  if (what == "coef") tex_coef(r$est, r$p) else tex_se(r$se)
}
for (g in groups) {
  rows <- list()
  for (y in names(raw_outcomes)) {
    rows[[length(rows) + 1]] <- c(raw_outcomes[[y]], sapply(names(transforms), \(tr) cell(g, y, tr, "coef")))
    rows[[length(rows) + 1]] <- c("", sapply(names(transforms), \(tr) cell(g, y, tr, "se")))
  }
  ctry <- filter(res, group == g)$n_countries
  rows <- c(rows, list("\\midrule", c("Countries", rep(paste0(min(ctry), "--", max(ctry)), 4))))
  write_tex_table(
    file    = file.path(out_dir, paste0("did_transformations_", group_file[[g]], ".tex")),
    caption = paste0("DiD under outcome transformations, adjusted: ", g, " vs.\\ controls"),
    label   = paste0("tab:crown_did_transform_", group_file[[g]]),
    head    = list(c("", "(1)", "(2)", "(3)", "(4)"),
                   c("", "Raw", "Winsorized 1/99", "log(x+1)", "IHS")),
    rows    = rows,
    notes   = paste(
      "Each cell is the DiD coefficient (treatment group $\\times$ Post) from Crown et al.\\ (2023)'s",
      "pooled OLS with the 11 baseline (2003) covariates, with the outcome transformed as in the",
      "column header. Raw: the outcome as reported (GDP per capita in constant 2022 US dollars, GDP in",
      "billions, growth and employment in \\%, disengagement in percentage points). Winsorized: values",
      "outside the 1st--99th percentiles of the full panel set to those percentiles. log(x+1): $\\ln(1 + y)$;",
      "for growth rates, $100 \\ln(1 + y/100)$, the growth rate in log points. IHS:",
      "inverse hyperbolic sine, $\\sinh^{-1}(y) = \\ln(y + \\sqrt{y^2 + 1})$, defined for negative values.",
      "With log(x+1) or IHS, coefficients on levels are approximately proportional changes. Sample 1990--2022; controls",
      "are the 64 comparison countries valid through 2022. Standard errors clustered by country.",
      "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lcccc", resize = FALSE)
}

# ---- LaTeX: parallel-trends test for each transformation, one table per group ----
fmt_pre <- function(p) {
  s <- if (p < 0.001) "$<$0.001" else sprintf("%.3f", p)
  if (p < 0.05) paste0("\\textbf{", s, "}") else s
}
for (g in groups) {
  rows <- list()
  for (panel_name in c("unadj", "adj")) {
    title <- if (panel_name == "unadj") "Panel A: unadjusted (country and year fixed effects)"
             else "Panel B: adjusted (plus the 11 baseline covariates $\\times$ year)"
    if (panel_name == "adj") rows[[length(rows) + 1]] <- "\\midrule"
    rows[[length(rows) + 1]] <- sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", title)
    for (y in names(raw_outcomes)) {
      r <- filter(res, group == g, outcome == y)
      rows[[length(rows) + 1]] <- c(raw_outcomes[[y]], sapply(names(transforms), \(tr)
        fmt_pre(r[[paste0("pre_p_", panel_name)]][r$transform == tr])))
    }
  }
  write_tex_table(
    file    = file.path(out_dir, paste0("pretrends_transformations_", group_file[[g]], ".tex")),
    caption = paste0("Parallel-trends test under outcome transformations: ", g, " vs.\\ controls"),
    label   = paste0("tab:crown_pretrends_transform_", group_file[[g]]),
    head    = list(c("", "(1)", "(2)", "(3)", "(4)"),
                   c("", "Raw", "Winsorized 1/99", "log(x+1)", "IHS")),
    rows    = rows,
    notes   = paste(
      "Each cell is the $p$-value of a joint F-test that all 13 pre-2004 event-study coefficients",
      "(1990--2002, relative to 2003) are zero, with the outcome transformed as in the column header",
      paste0("(transformations as in Table~\\ref{tab:crown_did_transform_", group_file[[g]], "})."),
      "Bold: $p < 0.05$, parallel trends rejected at 5\\%.",
      "Event studies with country and year fixed effects (Panel A), plus each baseline covariate",
      "interacted with year and income group $\\times$ year fixed effects (Panel B).",
      "Sample 1990--2022; controls are the 64 comparison countries valid through 2022.",
      "Standard errors clustered by country."),
    align   = "lcccc", resize = FALSE)
}

cat("Transformations done. Output in", out_dir, "\n")

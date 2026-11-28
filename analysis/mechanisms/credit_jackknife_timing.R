# Credit channel: jackknife, timing, and fiscal outcomes
# Author: Anthony Babu
#
# Three follow-ups on the savings / investment / FDI results (credit_baseyear_triple.R,
# savings_panel_triple.R), all with the triple interaction or the PEPFAR dummy from
# econ_setup.R:
#
# 1. Jackknife. Several credit-channel coefficients looked driven by Equatorial Guinea.
#    Every triple-interaction regression for FDI, gross domestic savings and gross fixed
#    capital formation (as shares of current GDP and of base-year GDP; raw, winsorized and
#    winsorized IHS; All PEPFAR and COP; unadjusted and adjusted) is re-estimated once for
#    each country, leaving it out. Reported: the share of runs significant at 10% and the
#    countries whose removal turns a significant result insignificant (or the reverse).
#
# 2. Timing. If longer lives raise saving and investment only gradually (through age
#    structure, accumulated assets), the effect should appear late. Post is split into
#    three periods, 2004-09, 2010-15 and 2016-22, each interacted with the PEPFAR dummy:
#      Y_it = a_i + g_t + sum_p b_p (PEPFAR_i x Period_p) [+ X_i x Post_t] + e_it
#    for the same outcomes (winsorized), with and without Equatorial Guinea.
#
# 3. Fiscal and financial outcomes. Health aid can work through government budgets and the
#    balance of payments rather than private saving. Same PEPFAR-dummy and triple-
#    interaction specifications for domestic government health spending, external share of
#    health spending, tax revenue, government and household consumption, private credit,
#    the current account, net ODA and external debt (WDI, data_pull/pull_wdi_mechanisms.R),
#    winsorized, without Equatorial Guinea. The two health-spending series start in 2000.
#
# Run from the repository root:  Rscript analysis/mechanisms/credit_jackknife_timing.R
# Writes to tables/crown_panel/mechanisms/credit_followups/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/mechanisms/credit_followups"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mech    <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)
panel_m <- add_baseyear_credit(panel, mech)

credit_outcomes <- c(fdi_pct_gdp     = "FDI (\\% of GDP)",
                     savings_pct_gdp = "Gross domestic savings (\\% of GDP)",
                     gfcf_pct_gdp    = "Fixed capital formation (\\% of GDP)",
                     fdi_gdp0        = "FDI (\\% of base-year GDP)",
                     dsav_gdp0       = "Gross domestic savings (\\% of base-year GDP)",
                     gfcf_gdp0       = "Fixed capital formation (\\% of base-year GDP)")
fiscal_outcomes <- c(dom_hexp_gdp    = "Domestic govt.\\ health spending (\\% of GDP)",
                     ext_hexp_che    = "External share of health spending (\\%)",
                     tax_gdp         = "Tax revenue (\\% of GDP)",
                     govcons_gdp     = "Government consumption (\\% of GDP)",
                     hhcons_gdp      = "Household consumption (\\% of GDP)",
                     credit_priv_gdp = "Private credit (\\% of GDP)",
                     cab_gdp         = "Current account balance (\\% of GDP)",
                     oda_gni         = "Net ODA received (\\% of GNI)",
                     ext_debt_gni    = "External debt (\\% of GNI)")
transforms <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs      <- c(unadj = "Unadj.", adj = "Adj.")

for (y in c(names(credit_outcomes), names(fiscal_outcomes))) {
  w <- winsor(panel_m[[y]])
  panel_m[[paste0(y, "__raw")]]    <- panel_m[[y]]
  panel_m[[paste0(y, "__win")]]    <- w
  panel_m[[paste0(y, "__winihs")]] <- asinh(w)
}

triple_coef <- function(d, yvar, adjusted) {
  m <- run_twfe_triple_post(d, yvar, adjusted = adjusted)
  list(m = m, ct = coeftable(m)["hiv_treat_post", ])
}
fmt_list <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"

# =====================================================
# 1. Jackknife of the triple interaction
# =====================================================
jk_runs <- list(); jk_full <- list()
for (g in groups) {
  d0 <- make_sample(g, data = panel_m) |> filter(!is.na(bl_hiv_prev))
  for (y in names(credit_outcomes)) for (tr in names(transforms)) for (sp in names(specs)) {
    yvar <- paste0(y, "__", tr)
    f0 <- triple_coef(d0, yvar, sp == "adj")
    key <- tibble(outcome = y, transform = tr, group = g, spec = sp)
    jk_full[[length(jk_full) + 1]] <- bind_cols(key, tibble(est = f0$ct[["Estimate"]], p = f0$ct[["Pr(>|t|)"]]))
    for (ctry in sort(unique(d0$iso3c[obs(f0$m)]))) {
      f <- triple_coef(filter(d0, iso3c != ctry), yvar, sp == "adj")
      jk_runs[[length(jk_runs) + 1]] <- bind_cols(key, tibble(dropped = ctry, est = f$ct[["Estimate"]],
                                                              p = f$ct[["Pr(>|t|)"]]))
    }
  }
}
jk_runs <- bind_rows(jk_runs)
jk_summ <- jk_runs |>
  left_join(bind_rows(jk_full), by = c("outcome", "transform", "group", "spec"), suffix = c("", "_full")) |>
  group_by(outcome, transform, group, spec) |>
  summarise(est_full = first(est_full), p_full = first(p_full), n_runs = n(),
            share_sig10 = mean(p < 0.10),
            flips = fmt_list(dropped[(p < 0.10) != (first(p_full) < 0.10)]),
            est_no_gnq = est[dropped == "GNQ"][1], p_no_gnq = p[dropped == "GNQ"][1], .groups = "drop")
write_csv(jk_runs, file.path(out_dir, "credit_jackknife_runs.csv"))
write_csv(jk_summ, file.path(out_dir, "credit_jackknife_summary.csv"))

# LaTeX: raw and winsorized IHS for each outcome, All PEPFAR and COP, unadjusted and adjusted
pct <- function(x) sprintf("%.0f\\%%", 100 * x)
cols <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
rows <- list()
for (y in names(credit_outcomes)) for (tr in c("raw", "winihs")) {
  r <- lapply(seq_len(nrow(cols)), \(j) filter(jk_summ, outcome == y, transform == tr,
                                               group == cols$group[j], spec == cols$spec[j]))
  rows[[length(rows) + 1]] <- c(sprintf("%s, %s", credit_outcomes[[y]], tolower(transforms[[tr]])),
                                sapply(r, \(x) sprintf("%s / %s", tex_coef(x$est_full, x$p_full),
                                                       pct(x$share_sig10))))
  rows[[length(rows) + 1]] <- c("\\quad without Equatorial Guinea",
                                sapply(r, \(x) tex_coef(x$est_no_gnq, x$p_no_gnq)))
}
write_tex_table(
  file    = file.path(out_dir, "credit_jackknife.tex"),
  caption = "Credit channel: jackknife of the HIV $\\times$ PEPFAR $\\times$ Post coefficient",
  label   = "tab:crown_credit_jackknife",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR sample}", "\\multicolumn{2}{c}{COP sample}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                 c("Outcome", rep(c("Unadj.", "Adj."), 2))),
  rows    = rows,
  notes   = paste(
    "First row of each pair: full-sample triple-interaction coefficient (HIV prevalence 2003 $\\times$",
    "PEPFAR $\\times$ Post, as in Section 5.3), followed by the share of leave-one-country-out",
    "regressions in which it is significant at 10\\%. Second row: the coefficient when Equatorial",
    "Guinea alone is left out. Each regression is re-estimated once for every country in its sample.",
    "Raw and winsorized IHS outcomes shown; the winsorized results are in the replication files.",
    "Country and year fixed effects; adjusted columns add the other 10 baseline covariates",
    "$\\times$ Post. Standard errors clustered by country. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align = "lcccc", resize = TRUE)

# =====================================================
# 2. Timing: PEPFAR effect by period
# =====================================================
periods <- c(p1 = "2004--09", p2 = "2010--15", p3 = "2016--22")
panel_m <- panel_m |>
  mutate(p1 = as.integer(year >= 2004 & year <= 2009),
         p2 = as.integer(year >= 2010 & year <= 2015),
         p3 = as.integer(year >= 2016))
run_periods <- function(d, yvar, adjusted) {
  d <- mutate(d, tp1 = treat * p1, tp2 = treat * p2, tp3 = treat * p3)
  rhs <- c("tp1", "tp2", "tp3")
  if (adjusted) rhs <- c(rhs, paste0(bl_cont, ":post"), "i(bl_wb_income, post, ref = 'L')")
  feols(as.formula(paste(yvar, "~", paste(rhs, collapse = " + "), "| iso3c + year")),
        data = d, cluster = ~iso3c)
}
tm <- list()
for (s in c("full", "nognq")) for (g in groups) {
  d <- make_sample(g, data = panel_m)
  if (s == "nognq") d <- filter(d, iso3c != "GNQ")
  for (y in names(credit_outcomes)) for (sp in names(specs)) {
    ct <- coeftable(run_periods(d, paste0(y, "__win"), sp == "adj"))
    for (k in seq_along(periods)) {
      v <- paste0("tp", k)
      tm[[length(tm) + 1]] <- tibble(sample = s, group = g, outcome = y, spec = sp, period = names(periods)[k],
                                     est = ct[v, "Estimate"], se = ct[v, "Std. Error"], p = ct[v, "Pr(>|t|)"])
    }
  }
}
tm <- bind_rows(tm)
write_csv(tm, file.path(out_dir, "credit_timing.csv"))

rows <- list()
for (s in c("nognq", "full")) {
  if (s == "full") rows[[length(rows) + 1]] <- "\\midrule"
  rows[[length(rows) + 1]] <- sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\",
                                      if (s == "nognq") "Excluding Equatorial Guinea" else "All countries")
  for (y in names(credit_outcomes)) {
    r <- lapply(c("unadj", "adj"), \(sp) lapply(names(periods), \(pp)
      filter(tm, sample == s, group == "All PEPFAR", outcome == y, spec == sp, period == pp)))
    rows[[length(rows) + 1]] <- c(credit_outcomes[[y]], unlist(lapply(r, \(rr) sapply(rr, \(x) tex_coef(x$est, x$p)))))
    rows[[length(rows) + 1]] <- c("", unlist(lapply(r, \(rr) sapply(rr, \(x) tex_se(x$se)))))
  }
}
write_tex_table(
  file    = file.path(out_dir, "credit_timing.tex"),
  caption = "Credit channel over time: PEPFAR effect by period, All PEPFAR vs.\\ controls",
  label   = "tab:crown_credit_timing",
  head    = list(c("", "\\multicolumn{3}{c}{Unadjusted}", "\\multicolumn{3}{c}{Adjusted}"),
                 "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
                 c("Outcome", rep(periods, 2))),
  rows    = rows,
  notes   = paste(
    "Each row is one regression of the winsorized outcome on PEPFAR $\\times$ (2004--09),",
    "PEPFAR $\\times$ (2010--15) and PEPFAR $\\times$ (2016--22), with country and year fixed",
    "effects; each coefficient is the change relative to 1990--2003 in that period. Adjusted",
    "columns add the 11 baseline covariates $\\times$ Post. If longer life expectancy raises saving",
    "and investment gradually, the later coefficients should be larger. 90 PEPFAR recipients vs.\\",
    "64 comparison countries; standard errors clustered by country.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align = "lcccccc", resize = TRUE)

# =====================================================
# 3. Fiscal and financial outcomes
# =====================================================
fs <- list()
for (g in groups) {
  d <- make_sample(g, data = panel_m) |> filter(iso3c != "GNQ")
  for (y in names(fiscal_outcomes)) for (sp in names(specs)) {
    yvar <- paste0(y, "__win")
    m1 <- run_twfe_post(d, yvar, tvar = "treat", adjusted = sp == "adj")
    c1 <- coeftable(m1)["treat:post", ]
    m3 <- run_twfe_triple_post(filter(d, !is.na(bl_hiv_prev)), yvar, adjusted = sp == "adj")
    c3 <- coeftable(m3)["hiv_treat_post", ]
    fs[[length(fs) + 1]] <- tibble(
      outcome = y, group = g, spec = sp,
      did = c1[["Estimate"]], did_se = c1[["Std. Error"]], did_p = c1[["Pr(>|t|)"]],
      triple = c3[["Estimate"]], triple_se = c3[["Std. Error"]], triple_p = c3[["Pr(>|t|)"]],
      n_countries = n_countries(m1, d))
  }
}
fs <- bind_rows(fs)
write_csv(fs, file.path(out_dir, "fiscal_outcomes.csv"))

rows <- list()
for (y in names(fiscal_outcomes)) {
  r <- lapply(seq_len(nrow(cols)), \(j) filter(fs, outcome == y, group == cols$group[j], spec == cols$spec[j]))
  rows[[length(rows) + 1]] <- c(fiscal_outcomes[[y]], sapply(r, \(x) tex_coef(x$did, x$did_p)))
  rows[[length(rows) + 1]] <- c("", sapply(r, \(x) tex_se(x$did_se)))
}
write_tex_table(
  file    = file.path(out_dir, "fiscal_outcomes.tex"),
  caption = "Fiscal and financial outcomes after PEPFAR, two-way fixed effects",
  label   = "tab:crown_fiscal",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR}", "\\multicolumn{2}{c}{COP}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                 c("Outcome", rep(c("Unadj.", "Adj."), 2))),
  rows    = rows,
  notes   = paste(
    "Coefficient on PEPFAR $\\times$ Post. The triple-interaction versions are in the replication",
    "files but are not reported: with Equatorial Guinea excluded, comparison countries have almost no",
    "variation in HIV prevalence, so the comparison-group HIV gradient is poorly identified.",
    "Outcomes from the World Bank WDI, winsorized at the 1st and 99th percentiles; the two",
    "health-spending series start in 2000, so their pre-period is 2000--03. Equatorial Guinea is",
    "excluded throughout. Country and year fixed effects; adjusted columns add the baseline",
    "covariates $\\times$ Post. Standard errors clustered by country.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align = "lcccc", resize = FALSE, size = "\\small")

cat("Credit-channel follow-ups done:", nrow(jk_runs), "jackknife regressions. Output in", out_dir, "\n")

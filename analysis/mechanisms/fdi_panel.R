# Investment mechanism: did FDI, saving and investment rise after PEPFAR?
# Author: Anthony Babu
#
# Mechanism: PEPFAR lowers adult mortality -> people expect longer working lives -> more
# saving and investment, and foreign investors see a healthier workforce. If so, these
# should rise after 2004 in PEPFAR countries, and more so where HIV was worse.
#
# Same panel (TWFE) specification as the growth robustness check
# (analysis/replication/twfe_growth_robustness.R), with a mechanism variable as outcome:
#   Y_it = a_i + g_t + b (T_i x Post_t) [+ X_i x Post_t] + e_it,   Post = 1 from 2004
#   T_i = the 0/1 PEPFAR (or COP) dummy, or HIV prevalence 2003 (bl_hiv_prev)
#   Adjusted adds the 11 baseline (2003) covariates x Post (10 in the HIV version)
#
# Outcomes (pulled by data_pull/pull_wdi_mechanisms.R):
#   fdi_pct_gdp      FDI, net inflows (% of GDP)
#   fdi_pc           FDI, net inflows per person (constant 2022 US$)
#   savings_pct_gdp  Gross domestic savings (% of GDP): GDP minus all consumption, so it
#                    excludes aid and remittances (gross national savings would count
#                    PEPFAR itself as income)
#   gfcf_pct_gdp     Gross fixed capital formation (% of GDP): all investment in
#                    machines, buildings and infrastructure, domestic and foreign
# Both scalings of FDI are reported because if PEPFAR raised GDP, a % of GDP measure can
# stay flat even when FDI rises.
#
# Transformations: Raw; Winsorized (1st/99th percentile of the full panel); Winsorized,
# IHS. FDI is lumpy (one mining deal can exceed 30% of GDP) and negative in years of
# disinvestment, so logs are not used; IHS is defined for negative values.
#
# Run from the repository root:  Rscript analysis/mechanisms/fdi_panel.R
# Writes to tables/crown_panel/mechanisms/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/mechanisms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mech <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)

mech_outcomes <- c(fdi_pct_gdp     = "FDI, net inflows (\\% of GDP)",
                   fdi_pc          = "FDI, net inflows per person (constant 2022 US\\$)",
                   savings_pct_gdp = "Gross domestic savings (\\% of GDP)",
                   gfcf_pct_gdp    = "Gross fixed capital formation (\\% of GDP)")
transforms <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs      <- c(unadj = "Unadjusted", adj = "Adjusted")

panel_m <- left_join(panel, mech, by = c("iso3c", "year"))
for (y in names(mech_outcomes)) {
  w <- winsor(panel_m[[y]])
  panel_m[[paste0(y, "__raw")]]    <- panel_m[[y]]
  panel_m[[paste0(y, "__win")]]    <- w
  panel_m[[paste0(y, "__winihs")]] <- asinh(w)
}

res <- list()
for (g in groups) {
  d <- make_sample(g, data = panel_m)
  for (tv in names(treatments)) for (y in names(mech_outcomes))
    for (tr in names(transforms)) for (sp in names(specs)) {
      m  <- run_twfe_post(d, paste0(y, "__", tr), tvar = tv, adjusted = sp == "adj")
      ct <- coeftable(m)[paste0(tv, ":post"), ]
      res[[length(res) + 1]] <- tibble(
        outcome = y, group = g, treatment = tv, transform = tr, spec = sp,
        est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
        n_obs = nobs(m), n_countries = n_countries(m, d))
    }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "fdi_panel.csv"))

# ---- LaTeX: one table. Columns: treatment x group x spec. Panels: one per outcome ----
cols <- expand.grid(spec = names(specs), group = groups, treatment = names(treatments),
                    stringsAsFactors = FALSE)
pick <- function(y, tr, j) filter(res, outcome == y, transform == tr, treatment == cols$treatment[j],
                                  group == cols$group[j], spec == cols$spec[j])
rows <- list()
for (k in seq_along(mech_outcomes)) {
  y <- names(mech_outcomes)[k]
  if (k > 1) rows[[length(rows) + 1]] <- "\\midrule"
  rows[[length(rows) + 1]] <- sprintf("\\multicolumn{9}{l}{\\textit{Panel %s: %s}} \\\\",
                                      LETTERS[k], mech_outcomes[[y]])
  for (tr in names(transforms)) {
    r <- lapply(seq_len(nrow(cols)), \(j) pick(y, tr, j))
    rows[[length(rows) + 1]] <- c(transforms[[tr]], sapply(r, \(x) tex_coef(x$est, x$p)))
    rows[[length(rows) + 1]] <- c("", sapply(r, \(x) tex_se(x$se)))
  }
  n <- sapply(seq_len(nrow(cols)), \(j) {
    v <- filter(res, outcome == y, treatment == cols$treatment[j], group == cols$group[j],
                spec == cols$spec[j])$n_countries
    if (min(v) == max(v)) as.character(v[1]) else paste0(min(v), "--", max(v))
  })
  rows[[length(rows) + 1]] <- c("Countries", n)
}

write_tex_table(
  file    = file.path(out_dir, "fdi_panel.tex"),
  caption = "Investment mechanism: FDI, saving and investment after PEPFAR, two-way fixed effects",
  label   = "tab:crown_mech_fdi",
  head    = list(c("", sprintf("\\multicolumn{4}{c}{Treatment: %s}", treatments)),
                 "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}",
                 c("", rep(c("\\multicolumn{2}{c}{All PEPFAR sample}",
                             "\\multicolumn{2}{c}{COP sample}"), 2)),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
                 c("", rep(c("Unadj.", "Adj."), 4)),
                 c("", sprintf("(%d)", 1:8))),
  rows    = rows,
  notes   = paste(
    "Each cell is the coefficient on treatment $\\times$ Post from a separate two-way fixed effects",
    "regression with country and year fixed effects; Post = 1 from 2004 for every country.",
    "Columns (1)--(4): treatment is the 0/1 PEPFAR (or COP) group dummy. Columns (5)--(8): treatment",
    "is HIV prevalence among ages 15--49 in 2003 (percentage points), so the coefficient is the",
    "post-2004 change per percentage point of baseline HIV prevalence.",
    "Adjusted columns add the 11 baseline (2003) covariates of Crown et al.\\ (2023), each interacted",
    "with Post (columns (6) and (8) leave out HIV prevalence, which is the treatment there).",
    "Outcomes from the World Bank WDI; FDI per person is net inflows in current US\\$ deflated to 2022",
    "dollars with the US GDP deflator and divided by UN WPP population. Winsorized: values outside the",
    "1st--99th percentiles of the full panel set to those percentiles. IHS: $\\sinh^{-1}(w)$ of the",
    "winsorized value, defined for the negative values that FDI and savings take in some years.",
    "Sample 1990--2022; controls are the 64 comparison countries valid through 2022.",
    "Standard errors clustered by country in parentheses.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccccccc", resize = TRUE)

cat("FDI / investment mechanism panel done. Output in", out_dir, "\n")

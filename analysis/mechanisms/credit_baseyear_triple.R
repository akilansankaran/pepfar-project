# Credit channel with fixed base-year denominators
# Author: Anthony Babu
#
# The mechanism tables divide each outcome by the same year's GDP or population. PEPFAR
# changes both: lower mortality raises population, and GDP grew in PEPFAR countries after
# 2004. A share of GDP or a per-person figure can therefore move even if saving and
# investment do not. Here every outcome is divided by the country's GDP or population in
# a fixed base year instead, the same number in every year, so changes in the measure come
# only from changes in the numerator.
#
# Base year: the first year of the panel, 1990, the earliest pre-PEPFAR year in the data.
# Countries whose GDP starts later (e.g. Eritrea 1992, the Baltic states and Serbia
# 1993-95) use their first year with data, as long as it is before 2004; countries whose
# GDP series starts in 2004 or later (Kosovo) are dropped, since their base would already
# be affected by PEPFAR. Population is available from 1990 for every country.
#
#   Share of base-year GDP (constant 2022 US$ in numerator and denominator):
#     fdi_gdp0      FDI net inflows        = fdi_pct_gdp     x GDP_t / GDP_base
#     dsav_gdp0     Gross domestic savings = savings_pct_gdp x GDP_t / GDP_base
#     gfcf_gdp0     Fixed capital formation= gfcf_pct_gdp    x GDP_t / GDP_base
#   Per person of the base year (constant 2022 US$):
#     fdi_pop0      FDI net inflows        = fdi_pc       x pop_t / pop_base
#     gsav_pop0     Gross savings          = gross_sav_pc x pop_t / pop_base
#     gfcf_pop0     Fixed capital formation= gfcf_pc      x pop_t / pop_base
# The WDI shares of GDP (data_pull/pull_wdi_mechanisms.R) times the panel's GDP in
# constant 2022 US$ give the numerator in constant 2022 US$, because both are deflated
# with the same US GDP deflator.
#
# Specification: the triple interaction used for the credit channel,
#   Y_it = a_i + g_t + b3 (HIV_i x PEPFAR_i x Post_t) + b1 (HIV_i x Post_t)
#                    + b2 (PEPFAR_i x Post_t) [+ X_i x Post_t] + e_it
# (run_twfe_triple_post in econ_setup.R); raw, winsorized and winsorized IHS outcomes;
# all countries and without Equatorial Guinea.
#
# Run from the repository root:  Rscript analysis/mechanisms/credit_baseyear_triple.R
# Writes to tables/crown_panel/mechanisms/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/mechanisms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mech <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)

# Each country's base-year GDP and population, and the six outcomes (econ_setup.R)
panel_m <- add_baseyear_credit(panel, mech)
write_csv(attr(panel_m, "bases"), file.path(out_dir, "credit_baseyear_bases.csv"))

tables <- list(
  gdp0 = list(file = "credit_baseyear_gdp_triple", label = "tab:crown_mech_gdp0_triple",
              caption = "Credit channel as a share of base-year GDP",
              unit = "percent of the country's GDP in the base year (1990 or first year with data)",
              outcomes = c(fdi_gdp0  = "FDI, net inflows (\\% of base-year GDP)",
                           dsav_gdp0 = "Gross domestic savings (\\% of base-year GDP)",
                           gfcf_gdp0 = "Gross fixed capital formation (\\% of base-year GDP)")),
  pop0 = list(file = "credit_baseyear_pop_triple", label = "tab:crown_mech_pop0_triple",
              caption = "Credit channel per person of the base-year population",
              unit = "constant 2022 US\\$ per person of the country's 1990 population",
              outcomes = c(fdi_pop0  = "FDI, net inflows (US\\$ per base-year person)",
                           gsav_pop0 = "Gross savings (US\\$ per base-year person)",
                           gfcf_pop0 = "Gross fixed capital formation (US\\$ per base-year person)")))
all_outcomes <- unlist(lapply(tables, `[[`, "outcomes"))
names(all_outcomes) <- sub("^[a-z0-9]+\\.", "", names(all_outcomes))

transforms   <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs        <- c(unadj = "Unadjusted", adj = "Adjusted")
samples_list <- c(full = "All countries", nognq = "Excluding Equatorial Guinea")

for (y in names(all_outcomes)) {
  w <- winsor(panel_m[[y]])
  panel_m[[paste0(y, "__raw")]]    <- panel_m[[y]]
  panel_m[[paste0(y, "__win")]]    <- w
  panel_m[[paste0(y, "__winihs")]] <- asinh(w)
}

res <- list()
for (s in names(samples_list)) for (g in groups) {
  d <- make_sample(g, data = panel_m) |> filter(!is.na(bl_hiv_prev))
  if (s == "nognq") d <- filter(d, iso3c != "GNQ")
  for (y in names(all_outcomes)) for (tr in names(transforms)) for (sp in names(specs)) {
    m  <- run_twfe_triple_post(d, paste0(y, "__", tr), adjusted = sp == "adj")
    ct <- coeftable(m)["hiv_treat_post", ]
    res[[length(res) + 1]] <- tibble(
      outcome = y, sample = s, group = g, transform = tr, spec = sp,
      est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
      n_obs = nobs(m), n_countries = n_countries(m, d), !!!as.list(triple_parts(m)))
  }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "credit_baseyear_triple.csv"))

# ---- LaTeX: one table per denominator. Columns: group x spec. Panels: outcome x sample ----
cols <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
pick <- function(y, s, tr, j) filter(res, outcome == y, sample == s, transform == tr,
                                     group == cols$group[j], spec == cols$spec[j])
for (tb in tables) {
  rows <- list()
  for (k in seq_along(tb$outcomes)) {
    y <- names(tb$outcomes)[k]
    if (k > 1) rows[[length(rows) + 1]] <- "\\midrule"
    rows[[length(rows) + 1]] <- sprintf("\\multicolumn{5}{l}{\\textbf{Panel %s: %s}} \\\\",
                                        LETTERS[k], tb$outcomes[[y]])
    for (s in names(samples_list)) {
      rows[[length(rows) + 1]] <- sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", samples_list[[s]])
      for (tr in names(transforms)) {
        r <- lapply(seq_len(nrow(cols)), \(j) pick(y, s, tr, j))
        rows[[length(rows) + 1]] <- c(transforms[[tr]], sapply(r, \(x) tex_coef(x$est, x$p)))
        rows[[length(rows) + 1]] <- c("", sapply(r, \(x) tex_se(x$se)))
      }
      n <- sapply(seq_len(nrow(cols)), \(j) {
        v <- filter(res, outcome == y, sample == s, group == cols$group[j], spec == cols$spec[j])$n_countries
        if (min(v) == max(v)) as.character(v[1]) else paste0(min(v), "--", max(v))
      })
      rows[[length(rows) + 1]] <- c("Countries", n)
    }
  }
  write_tex_table(
    file    = file.path(out_dir, paste0(tb$file, ".tex")),
    caption = paste0(tb$caption, ": HIV prevalence (2003) $\\times$ PEPFAR $\\times$ Post"),
    label   = tb$label,
    head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR sample}", "\\multicolumn{2}{c}{COP sample}"),
                   "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                   c("", rep(c("Unadj.", "Adj."), 2)),
                   c("", sprintf("(%d)", 1:4))),
    rows    = rows,
    notes   = paste(
      "Each cell is the coefficient on HIV prevalence $\\times$ PEPFAR $\\times$ Post from a two-way",
      "fixed effects regression with country and year fixed effects that also includes HIV",
      "prevalence $\\times$ Post and PEPFAR $\\times$ Post; Post = 1 from 2004. HIV prevalence is",
      "among ages 15--49 in 2003; PEPFAR is the 0/1 group dummy (All PEPFAR or COP, by column).",
      "Adjusted columns add the other 10 baseline (2003) covariates of Crown et al.\\ (2023), each",
      "interacted with Post. Outcomes are in", tb$unit, "so the denominator is fixed for each",
      "country and cannot respond to PEPFAR. Base year: 1990, or the first later year with data",
      "if before 2004 (countries whose series starts in 2004 or later are dropped). Numerators",
      "from the World Bank WDI, in constant 2022 US\\$ (US GDP deflator). Winsorized: values",
      "outside the 1st--99th percentiles of the full panel set to those percentiles. IHS:",
      "$\\sinh^{-1}(w)$ of the winsorized value. The second sub-panel of each outcome drops",
      "Equatorial Guinea, the only control with 2003 HIV prevalence above 1\\% (3.8\\%).",
      "Sample 1990--2022; controls are the 64 comparison countries valid through 2022. Standard",
      "errors clustered by country in parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lcccc", resize = FALSE, size = "\\small")
}

cat("Credit channel with base-year denominators done. Output in", out_dir, "\n")

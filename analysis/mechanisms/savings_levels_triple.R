# Credit channel, not normed by GDP: gross savings and fixed capital formation per person
# Author: Anthony Babu
#
# Same triple-interaction TWFE as the other mechanism tables (savings_panel_triple.R):
#   Y_it = a_i + g_t + b3 (HIV_i x PEPFAR_i x Post_t) + b1 (HIV_i x Post_t)
#                    + b2 (PEPFAR_i x Post_t) [+ X_i x Post_t] + e_it
# (run_twfe_triple_post in econ_setup.R), with outcomes measured per person rather than as
# a share of GDP. GDP itself rose in PEPFAR countries after 2004, so a share of GDP can stay
# flat while saving or investment rises; per person in constant dollars avoids that.
#
#   gross_sav_pc  Gross savings per person, constant 2022 US$ (WDI NY.GNS.ICTR.CD / population)
#   gfcf_pc       Gross fixed capital formation per person, constant 2022 US$
#                 (WDI NE.GDI.FTOT.CD / population): all investment in machines, buildings
#                 and infrastructure, the "capital formation" measure in levels
# Both are pulled in current US$ (not local currency) so every country is in the same
# dollars, then deflated with the US GDP deflator like the rest of the panel
# (data_pull/pull_wdi_mechanisms.R).
#
# Transformations: Raw; Winsorized (1st/99th percentile of the full panel); Winsorized,
# IHS, which is close to a log for values in the hundreds of dollars, so the coefficient
# is roughly a proportional change, and is defined for the negative savings some
# countries report.
#
# Samples: all countries, and without Equatorial Guinea, the only control with 2003 HIV
# prevalence above 1% (3.8%). The controls' HIV gradient b1, against which b3 is
# measured, rests largely on it.
#
# Run from the repository root:  Rscript analysis/mechanisms/savings_levels_triple.R
# Writes to tables/crown_panel/mechanisms/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/mechanisms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mech <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)

level_outcomes <- c(gross_sav_pc = "Gross savings per person (constant 2022 US\\$)",
                    gfcf_pc      = "Gross fixed capital formation per person (constant 2022 US\\$)")
transforms   <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs        <- c(unadj = "Unadjusted", adj = "Adjusted")
samples_list <- c(full = "All countries", nognq = "Excluding Equatorial Guinea")

panel_m <- left_join(panel, mech, by = c("iso3c", "year"))
for (y in names(level_outcomes)) {
  w <- winsor(panel_m[[y]])
  panel_m[[paste0(y, "__raw")]]    <- panel_m[[y]]
  panel_m[[paste0(y, "__win")]]    <- w
  panel_m[[paste0(y, "__winihs")]] <- asinh(w)
}

res <- list()
for (s in names(samples_list)) for (g in groups) {
  d <- make_sample(g, data = panel_m) |> filter(!is.na(bl_hiv_prev))
  if (s == "nognq") d <- filter(d, iso3c != "GNQ")
  for (y in names(level_outcomes)) for (tr in names(transforms)) for (sp in names(specs)) {
    m  <- run_twfe_triple_post(d, paste0(y, "__", tr), adjusted = sp == "adj")
    ct <- coeftable(m)["hiv_treat_post", ]
    res[[length(res) + 1]] <- tibble(
      outcome = y, sample = s, group = g, transform = tr, spec = sp,
      est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
      n_obs = nobs(m), n_countries = n_countries(m, d), !!!as.list(triple_parts(m)))
  }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "savings_levels_triple.csv"))

# ---- LaTeX: one table. Columns: group x spec. Panels: outcome; sub-panels: sample ----
cols <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
pick <- function(y, s, tr, j) filter(res, outcome == y, sample == s, transform == tr,
                                     group == cols$group[j], spec == cols$spec[j])
rows <- list()
for (k in seq_along(level_outcomes)) {
  y <- names(level_outcomes)[k]
  if (k > 1) rows[[length(rows) + 1]] <- "\\midrule"
  rows[[length(rows) + 1]] <- sprintf("\\multicolumn{5}{l}{\\textbf{Panel %s: %s}} \\\\",
                                      LETTERS[k], level_outcomes[[y]])
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
  file    = file.path(out_dir, "savings_levels_triple.tex"),
  caption = paste("Credit channel in levels: gross savings and fixed capital formation per person,",
                  "HIV prevalence (2003) $\\times$ PEPFAR $\\times$ Post"),
  label   = "tab:crown_mech_levels_triple",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR sample}", "\\multicolumn{2}{c}{COP sample}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                 c("", rep(c("Unadj.", "Adj."), 2)),
                 c("", sprintf("(%d)", 1:4))),
  rows    = rows,
  notes   = paste(
    "Each cell is the coefficient on HIV prevalence $\\times$ PEPFAR $\\times$ Post from a two-way",
    "fixed effects regression with country and year fixed effects that also includes HIV",
    "prevalence $\\times$ Post and PEPFAR $\\times$ Post; Post = 1 from 2004. HIV prevalence is",
    "among ages 15--49 in 2003 (percentage points); PEPFAR is the 0/1 group dummy (All PEPFAR or",
    "COP, by column). The coefficient is how much more the post-2004 change in the outcome rises",
    "per percentage point of HIV prevalence in the treatment group than in the controls. Adjusted",
    "columns add the other 10 baseline (2003) covariates of Crown et al.\\ (2023), each interacted",
    "with Post. Outcomes from the World Bank WDI (gross savings NY.GNS.ICTR.CD, gross fixed",
    "capital formation NE.GDI.FTOT.CD), in current US\\$ deflated to 2022 dollars with the US GDP",
    "deflator and divided by UN WPP population. Winsorized: values outside the 1st--99th",
    "percentiles of the full panel set to those percentiles. IHS: $\\sinh^{-1}(w)$ of the",
    "winsorized value, roughly a proportional change. The second sub-panel drops Equatorial",
    "Guinea, the only control with 2003 HIV prevalence above 1\\% (3.8\\%); it has no gross savings",
    "data, so the two sub-panels of Panel A are identical. Sample 1990--2022;",
    "controls are the 64 comparison countries valid through 2022. Standard errors clustered by",
    "country in parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccc", resize = FALSE, size = "\\small")

cat("Savings / investment levels triple interaction done. Output in", out_dir, "\n")

# How robust is the TWFE growth effect? Winsorizing, logs, IHS, and dropping Asia
# Author: Anthony Babu
#
# The panel (TWFE) regressions find a significant, positive effect of PEPFAR on GDP
# per capita growth, unadjusted and adjusted. This script checks how sensitive that
# result is, for GDP per capita growth and GDP growth:
#
#   Specification (two-way fixed effects, country and year FE, Post = 1 from 2004):
#     Unadjusted  Y_it = a_i + g_t + b (PEPFAR_i x Post_t) + e_it
#     Adjusted    adds the 11 baseline (2003) covariates of Crown et al. (2023), each
#                 interacted with Post: X_i x Post_t. They are constant within a country,
#                 so on their own the country FE would absorb them; x Post lets each
#                 shift growth after 2004. Same form as the Kates panel
#                 (mortality_panel_twfe.R).
#
#   Outcome transformations (cut-offs from the full 1990-2022 panel):
#     Raw              growth in % as reported ("levels")
#     Winsorized       values below the 1st / above the 99th percentile set to those
#                      percentiles, which blunts outliers such as Equatorial Guinea's
#                      +140% growth in 1997
#     Winsorized, log  100 * log(1 + w/100) of the winsorized growth rate w: the growth
#                      rate in log points. log(1 + w) with w in % would drop every year
#                      with growth at or below -1%
#     Winsorized, IHS  asinh(w) = log(w + sqrt(w^2 + 1)) of the winsorized growth rate;
#                      like a log for large values, defined for zero and negatives
#
#   Treatment: the 0/1 PEPFAR (or COP) dummy x Post, and the triple interaction
#   HIV prevalence 2003 x PEPFAR x Post (run_twfe_triple_post in econ_setup.R), which asks
#   whether the post-2004 change rises with baseline HIV more in PEPFAR countries than in
#   the controls. HIV prevalence is adult (15-49) prevalence in 2003 (bl_hiv_prev, already %
#   of population), fixed before PEPFAR so PEPFAR's own effect on prevalence cannot feed
#   back into it; the same variable as the Kates trend breaks. PEPFAR money went mostly to
#   high-HIV countries, so if PEPFAR raised growth the effect should scale with it.
#
#   Samples: all PEPFAR (90) or COP (31) vs. the 64 comparison countries; and the same
#   samples excluding Asia, so that the 1997-98 Asian financial crisis, which hit just
#   before PEPFAR began, cannot drive the before/after comparison; and the same samples
#   excluding Equatorial Guinea, the only control with 2003 HIV prevalence above 1% (3.8%)
#   and an oil boom that ended just before 2004, which dominates the controls' HIV gradient
#   in the triple interaction.
#
# Run from the repository root:  Rscript analysis/replication/twfe_growth_robustness.R
# Writes to tables/crown_panel/growth_robustness/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/growth_robustness"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

growth_outcomes <- c(gdp_pc_growth = "GDP per capita growth",
                     gdp_growth    = "GDP growth")

# East, Southeast and South Asia (World Bank regions East Asia & Pacific and South Asia,
# without the Pacific islands). Central Asia and the Middle East stay in the sample.
# Treated: AFG BGD CHN IDN IND KHM LAO MMR NPL PAK PHL THA TLS VNM (COP: IDN IND KHM MMR VNM)
# Controls: BTN LKA MDV MNG MYS PRK
asia <- c("AFG", "BGD", "BTN", "CHN", "IDN", "IND", "KHM", "LAO", "LKA", "MDV",
          "MMR", "MNG", "MYS", "NPL", "PAK", "PHL", "PRK", "THA", "TLS", "VNM")

transforms <- c(raw = "Raw (levels)", win = "Winsorized",
                winlog = "Winsorized, log", winihs = "Winsorized, IHS")
panel_t <- panel
for (y in names(growth_outcomes)) {
  w <- winsor(panel[[y]])
  panel_t[[paste0(y, "__raw")]]    <- panel[[y]]
  panel_t[[paste0(y, "__win")]]    <- w
  panel_t[[paste0(y, "__winlog")]] <- 100 * log1p(w / 100)
  panel_t[[paste0(y, "__winihs")]] <- asinh(w)
}

samples_list <- c(full = "All countries", noasia = "Excluding Asia",
                  nognq = "Excluding Equatorial Guinea")
treatments_g <- c(treat  = "PEPFAR group $\\times$ Post",
                  triple = "HIV prevalence 2003 $\\times$ PEPFAR $\\times$ Post")
specs        <- c(unadj = "Unadjusted", adj = "Adjusted")

res <- list()
for (s in names(samples_list)) for (g in groups) {
  d <- make_sample(g, data = panel_t)
  if (s == "noasia") d <- filter(d, !iso3c %in% asia)
  if (s == "nognq")  d <- filter(d, iso3c != "GNQ")
  for (tv in names(treatments_g)) for (y in names(growth_outcomes))
    for (tr in names(transforms)) for (sp in names(specs)) {
      yv <- paste0(y, "__", tr)
      if (tv == "treat") {
        m  <- run_twfe_post(d, yv, tvar = "treat", adjusted = sp == "adj")
        ct <- coeftable(m)["treat:post", ]
      } else {
        m  <- run_twfe_triple_post(d, yv, adjusted = sp == "adj")
        ct <- coeftable(m)["hiv_treat_post", ]
      }
      parts <- if (tv == "triple") triple_parts(m) else
        c(b1 = NA, b1_se = NA, slope_treated = NA, slope_treated_se = NA)
      res[[length(res) + 1]] <- tibble(
        outcome = y, sample = s, group = g, treatment = tv, transform = tr, spec = sp,
        est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
        n_obs = nobs(m), n_countries = n_countries(m, d), !!!as.list(parts))
    }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "twfe_growth_robustness.csv"))

# ---- LaTeX: one table per outcome. Columns: treatment x group x spec. Panels: full / no Asia ----
cols <- expand.grid(spec = names(specs), group = groups, treatment = names(treatments_g),
                    stringsAsFactors = FALSE)
cell <- function(y, s, tr, what) {
  sapply(seq_len(nrow(cols)), \(j) {
    r <- filter(res, outcome == y, sample == s, transform == tr, treatment == cols$treatment[j],
                group == cols$group[j], spec == cols$spec[j])
    if (what == "coef") tex_coef(r$est, r$p) else tex_se(r$se)
  })
}
countries <- function(y, s) {
  sapply(seq_len(nrow(cols)), \(j) {
    n <- filter(res, outcome == y, sample == s, treatment == cols$treatment[j],
                group == cols$group[j], spec == cols$spec[j])$n_countries
    if (min(n) == max(n)) as.character(n[1]) else paste0(min(n), "--", max(n))
  })
}

for (y in names(growth_outcomes)) {
  rows <- list()
  for (s in names(samples_list)) {
    if (s != "full") rows[[length(rows) + 1]] <- "\\midrule"
    rows[[length(rows) + 1]] <- sprintf("\\multicolumn{9}{l}{\\textit{Panel %s: %s}} \\\\",
                                        LETTERS[match(s, names(samples_list))], samples_list[[s]])
    for (tr in names(transforms)) {
      rows[[length(rows) + 1]] <- c(transforms[[tr]], cell(y, s, tr, "coef"))
      rows[[length(rows) + 1]] <- c("", cell(y, s, tr, "se"))
    }
    rows[[length(rows) + 1]] <- c("Countries", countries(y, s))
  }
  write_tex_table(
    file    = file.path(out_dir, paste0("twfe_robustness_", y, ".tex")),
    caption = paste0("PEPFAR and ", growth_outcomes[[y]],
                     ": two-way fixed effects with the PEPFAR dummy and the HIV $\\times$ PEPFAR triple interaction, under outcome transformations"),
    label   = paste0("tab:crown_twfe_robust_", y),
    head    = list(c("", sprintf("\\multicolumn{4}{c}{%s}", treatments_g)),
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
      "Columns (1)--(4): coefficient on the 0/1 PEPFAR (or COP) group dummy $\\times$ Post.",
      "Columns (5)--(8): coefficient on HIV prevalence $\\times$ PEPFAR $\\times$ Post from a regression",
      "that also includes HIV prevalence $\\times$ Post and PEPFAR $\\times$ Post; it is how much more the",
      "post-2004 change in growth rises per percentage point of 2003 HIV prevalence (ages 15--49) in",
      "the treatment group than in the controls. Adjusted columns add the 11 baseline (2003)",
      "covariates of Crown et al.\\ (2023), each interacted with Post (columns (6) and (8) leave out",
      "HIV prevalence, which is already interacted with Post).",
      "Raw: growth in \\% as reported. Winsorized: values outside the 1st--99th percentiles of",
      "the full panel set to those percentiles. Winsorized, log: $100\\ln(1 + w/100)$ of the winsorized",
      "growth rate $w$, the growth rate in log points. Winsorized, IHS: $\\sinh^{-1}(w) = \\ln(w + \\sqrt{w^2+1})$.",
      "Panel B drops East, Southeast and South Asia (20 countries, 14 of them PEPFAR recipients and 5 COP)",
      "so that the 1997--98 Asian financial crisis cannot drive the before/after comparison.",
      "Panel C drops Equatorial Guinea, the only control with 2003 HIV prevalence above 1\\% (3.8\\%),",
      "whose oil boom ended just before 2004.",
      "Sample 1990--2022; controls are the 64 comparison countries valid through 2022.",
      "Standard errors clustered by country in parentheses.",
      "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lcccccccc", resize = TRUE)
}

cat("TWFE growth robustness done. Output in", out_dir, "\n")

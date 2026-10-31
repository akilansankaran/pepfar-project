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
#   Samples: all PEPFAR (90) or COP (31) vs. the 64 comparison countries; and the same
#   samples excluding Asia, so that the 1997-98 Asian financial crisis, which hit just
#   before PEPFAR began, cannot drive the before/after comparison.
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

winsor <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

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

samples_list <- c(full = "All countries", noasia = "Excluding Asia")
specs        <- c(unadj = "Unadjusted", adj = "Adjusted")

run_twfe <- function(d, yvar, adjusted) {
  rhs <- "treat:post"
  if (adjusted) {
    rhs <- c(rhs, paste0(bl_cont, ":post"), "i(bl_wb_income, post, ref = 'L')")
  }
  feols(as.formula(paste(yvar, "~", paste(rhs, collapse = " + "), "| iso3c + year")),
        data = d, cluster = ~iso3c)
}

res <- list()
for (s in names(samples_list)) for (g in groups) {
  d <- make_sample(g, data = panel_t)
  if (s == "noasia") d <- filter(d, !iso3c %in% asia)
  for (y in names(growth_outcomes)) for (tr in names(transforms)) for (sp in names(specs)) {
    yvar <- paste0(y, "__", tr)
    m <- run_twfe(d, yvar, adjusted = sp == "adj")
    ct <- coeftable(m)["treat:post", ]
    res[[length(res) + 1]] <- tibble(
      outcome = y, sample = s, group = g, transform = tr, spec = sp,
      est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
      n_obs = nobs(m), n_countries = n_countries(m, d),
      n_treated = length(unique(d$iso3c[obs(m)][d$treat[obs(m)] == 1])))
  }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "twfe_growth_robustness.csv"))

# ---- LaTeX: one table per outcome. Columns: group x spec. Panels: full / no Asia ----
cols <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
cell <- function(y, s, tr, what) {
  sapply(seq_len(nrow(cols)), \(j) {
    r <- filter(res, outcome == y, sample == s, transform == tr,
                group == cols$group[j], spec == cols$spec[j])
    if (what == "coef") tex_coef(r$est, r$p) else tex_se(r$se)
  })
}
countries <- function(y, s) {
  sapply(seq_len(nrow(cols)), \(j) {
    n <- filter(res, outcome == y, sample == s, group == cols$group[j], spec == cols$spec[j])$n_countries
    if (min(n) == max(n)) as.character(n[1]) else paste0(min(n), "--", max(n))
  })
}

for (y in names(growth_outcomes)) {
  rows <- list()
  for (s in names(samples_list)) {
    if (s == "noasia") rows[[length(rows) + 1]] <- "\\midrule"
    rows[[length(rows) + 1]] <- sprintf("\\multicolumn{5}{l}{\\textit{Panel %s: %s}} \\\\",
                                        if (s == "full") "A" else "B", samples_list[[s]])
    for (tr in names(transforms)) {
      rows[[length(rows) + 1]] <- c(transforms[[tr]], cell(y, s, tr, "coef"))
      rows[[length(rows) + 1]] <- c("", cell(y, s, tr, "se"))
    }
    rows[[length(rows) + 1]] <- c("Countries", countries(y, s))
  }
  write_tex_table(
    file    = file.path(out_dir, paste0("twfe_robustness_", y, ".tex")),
    caption = paste0("PEPFAR and ", growth_outcomes[[y]],
                     ": two-way fixed effects under outcome transformations and without Asia"),
    label   = paste0("tab:crown_twfe_robust_", y),
    head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR vs.\\ controls}",
                     "\\multicolumn{2}{c}{COP vs.\\ controls}"),
                   "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                   c("", "Unadjusted", "Adjusted", "Unadjusted", "Adjusted"),
                   c("", "(1)", "(2)", "(3)", "(4)")),
    rows    = rows,
    notes   = paste(
      "Each cell is the DiD coefficient (treatment group $\\times$ Post) from a separate two-way fixed",
      "effects regression with country and year fixed effects; Post = 1 from 2004 for every country.",
      "Adjusted columns add the 11 baseline (2003) covariates of Crown et al.\\ (2023), each interacted",
      "with Post. Raw: growth in \\% as reported. Winsorized: values outside the 1st--99th percentiles of",
      "the full panel set to those percentiles. Winsorized, log: $100\\ln(1 + w/100)$ of the winsorized",
      "growth rate $w$, the growth rate in log points. Winsorized, IHS: $\\sinh^{-1}(w) = \\ln(w + \\sqrt{w^2+1})$.",
      "Panel B drops East, Southeast and South Asia (20 countries, 14 of them PEPFAR recipients and 5 COP)",
      "so that the 1997--98 Asian financial crisis cannot drive the before/after comparison.",
      "Sample 1990--2022; controls are the 64 comparison countries valid through 2022.",
      "Standard errors clustered by country in parentheses.",
      "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lcccc", resize = FALSE)
}

cat("TWFE growth robustness done. Output in", out_dir, "\n")

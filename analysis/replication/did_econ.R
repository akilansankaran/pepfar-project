# Difference-in-differences for the Crown et al. (2023) economic outcomes,
# adjusted for all 11 baseline covariates, All PEPFAR and COP vs. controls, 1990-2022
# Author: Anthony Babu
#
# Two specifications for every outcome x treatment group:
#   Panel A (Crown's own model, pooled OLS, no fixed effects):
#     Y_it = b0 + b1 Post_t + b2 PEPFAR_i + b3 (PEPFAR_i x Post_t) + d'X_i + e_it
#   Panel B (two-way fixed effects, the same controls the adjusted event study uses):
#     Y_it = a_i + g_t + b (PEPFAR_i x Post_t) + X_i' d_t + e_it
#     a_i = country FE, g_t = year FE, X_i' d_t = each baseline covariate with its own
#     coefficient in every year (+ income group x year FE).
# b3 / b is the DiD estimate. SEs clustered by country.
#
# Run from the repository root:  Rscript analysis/replication/did_econ.R
# Writes LaTeX tables and a CSV to tables/crown_did/adjusted/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_did/adjusted"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. The two regressions
# =====================================================

run_did_crown <- function(d, y) {
  f <- as.formula(paste(y, "~ post + treat + post:treat +",
                        paste(bl_all, collapse = " + ")))
  feols(f, data = d, cluster = ~iso3c)
}

run_did_twfe <- function(d, y) {
  f <- as.formula(paste0(y, " ~ treat:post | iso3c + bl_wb_income^year + year[",
                         paste(bl_cont, collapse = ", "), "]"))
  feols(f, data = d, cluster = ~iso3c)
}

# One row per (model, term) with everything the tables need
grab <- function(m, d, y, g, spec, terms) {
  ct <- coeftable(m)
  tibble(outcome = y, group = g, spec = spec, term = terms,
         est = ct[terms, "Estimate"], se = ct[terms, "Std. Error"],
         p = ct[terms, "Pr(>|t|)"], n_obs = nobs(m),
         n_countries = n_countries(m, d), r2 = unname(r2(m, "r2")))
}

res <- list()
for (g in groups) for (y in names(outcomes)) {
  d <- samples[[g]]
  res[[length(res) + 1]] <- grab(run_did_crown(d, y), d, y, g, "Pooled OLS",
                                 c("post:treat", "treat", "post", "(Intercept)"))
  res[[length(res) + 1]] <- grab(run_did_twfe(d, y), d, y, g, "TWFE", "treat:post")
}
res <- bind_rows(res) |>
  mutate(term = recode(term, "post:treat" = "did", "treat:post" = "did"))
write_csv(res, file.path(out_dir, "did_adjusted.csv"))

# =====================================================
# 2. LaTeX tables
# =====================================================

ys <- names(outcomes)
cell <- function(g, spec, term, y, what) {
  r <- filter(res, group == g, spec == !!spec, term == !!term, outcome == y)
  switch(what,
         coef = tex_coef(r$est, r$p), se = tex_se(r$se),
         n = tex_n(r$n_obs), countries = tex_n(r$n_countries),
         r2 = sprintf("%.3f", r$r2))
}

did_notes <- paste(
  "Standard errors clustered by country in parentheses. Post $= 1$ for 2004",
  "onward for every country. Sample 1990--2022; controls are the 64 comparison",
  "countries valid through 2022",
  "(\\texttt{grp\\_control\\_kates}). All models include the 11 baseline (2003)",
  "covariates of Crown et al.\\ (2023) Table~1: GDP per capita, population, life",
  "expectancy, fertility, urban share, gross secondary enrollment, HIV prevalence,",
  "non-PEPFAR donor health aid per capita, domestic health expenditure per capita,",
  "US HIV aid before 2004, and World Bank income group. Pooled OLS enters them",
  "directly (no fixed effects). TWFE has country and year fixed",
  "effects and lets each covariate have its own coefficient in every year",
  "(income group enters as income group $\\times$ year fixed effects).",
  "Disengagement is in percentage points.",
  "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.")

# ---- Separate tables: one per treatment group, columns = the 8 outcomes ----
for (g in groups) {
  tl <- if (g == "COP") "COP" else "PEPFAR"
  panel_row <- function(txt) sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                                     length(ys) + 1, txt)
  rows <- list(
    panel_row("Panel A. Pooled OLS (Crown et al.\\ 2023 specification)"),
    c(paste(tl, "$\\times$ Post"), sapply(ys, \(y) cell(g, "Pooled OLS", "did", y, "coef"))),
    c("",                          sapply(ys, \(y) cell(g, "Pooled OLS", "did", y, "se"))),
    c(tl,                          sapply(ys, \(y) cell(g, "Pooled OLS", "treat", y, "coef"))),
    c("",                          sapply(ys, \(y) cell(g, "Pooled OLS", "treat", y, "se"))),
    c("Post",                      sapply(ys, \(y) cell(g, "Pooled OLS", "post", y, "coef"))),
    c("",                          sapply(ys, \(y) cell(g, "Pooled OLS", "post", y, "se"))),
    c("Observations",              sapply(ys, \(y) cell(g, "Pooled OLS", "did", y, "n"))),
    c("Countries",                 sapply(ys, \(y) cell(g, "Pooled OLS", "did", y, "countries"))),
    c("$R^2$",                     sapply(ys, \(y) cell(g, "Pooled OLS", "did", y, "r2"))),
    "\\midrule",
    panel_row("Panel B. Two-way fixed effects"),
    c(paste(tl, "$\\times$ Post"), sapply(ys, \(y) cell(g, "TWFE", "did", y, "coef"))),
    c("",                          sapply(ys, \(y) cell(g, "TWFE", "did", y, "se"))),
    c("Observations",              sapply(ys, \(y) cell(g, "TWFE", "did", y, "n"))),
    c("Countries",                 sapply(ys, \(y) cell(g, "TWFE", "did", y, "countries"))),
    c("$R^2$",                     sapply(ys, \(y) cell(g, "TWFE", "did", y, "r2")))
  )
  write_tex_table(
    file    = file.path(out_dir, paste0("did_adjusted_", group_file[[g]], ".tex")),
    caption = paste0("Difference-in-differences, adjusted: ", g, " vs.\\ controls"),
    label   = paste0("tab:did_", group_file[[g]]),
    head    = list(c("", sprintf("(%d)", seq_along(ys))), c("", outcomes_tex[ys])),
    rows    = rows,
    notes   = paste("Each column is a separate regression. Panel A: pooled OLS;",
                    "Panel B: two-way fixed effects (TWFE).", did_notes),
    align   = paste0("l", strrep("c", length(ys))))
}

# ---- Combined table: rows = outcomes, columns = All PEPFAR and COP side by side ----
combo_cols <- list(c("All PEPFAR", "Pooled OLS"), c("All PEPFAR", "TWFE"),
                   c("COP", "Pooled OLS"),        c("COP", "TWFE"))
rows <- list()
for (y in ys) {
  rows[[length(rows) + 1]] <- c(outcomes_tex[[y]],
                                sapply(combo_cols, \(k) cell(k[1], k[2], "did", y, "coef")))
  rows[[length(rows) + 1]] <- c("", sapply(combo_cols, \(k) cell(k[1], k[2], "did", y, "se")))
}
ctry_range <- sapply(combo_cols, \(k) {
  n <- filter(res, group == k[1], spec == k[2], term == "did")$n_countries
  paste0(min(n), "--", max(n))
})
rows <- c(rows, list("\\midrule", c("Countries (range across outcomes)", ctry_range)))
write_tex_table(
  file    = file.path(out_dir, "did_adjusted_combined.tex"),
  caption = "Difference-in-differences, adjusted: All PEPFAR and COP",
  label   = "tab:did_combined",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR vs.\\ controls}",
                   "\\multicolumn{2}{c}{COP vs.\\ controls}"),
                 "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
                 c("", "Pooled OLS", "TWFE", "Pooled OLS", "TWFE"),
                 c("", "(1)", "(2)", "(3)", "(4)")),
  rows    = rows,
  notes   = paste("Each cell is the DiD coefficient (treatment group $\\times$ Post)",
                  "from a separate regression; full results in",
                  "Tables~\\ref{tab:did_all_pepfar} and~\\ref{tab:did_cop}.", did_notes),
  align   = "lcccc", resize = FALSE)

cat("DiD done. Tables in", out_dir, "\n")

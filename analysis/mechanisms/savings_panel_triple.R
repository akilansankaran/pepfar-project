# Savings mechanism with a triple interaction: HIV prevalence x PEPFAR x Post
# Author: Rishad Vaghaiwalla
#
# Same regressions as Anthony's analysis/mechanisms/fdi_panel.R for gross domestic
# savings (% of GDP), with the same samples, transformations and controls, but the
# treatment term is HIV_i x PEPFAR_i x Post_t instead of HIV_i x Post_t:
#
#   Y_it = a_i + g_t + b3 (HIV_i x PEPFAR_i x Post_t) + b1 (HIV_i x Post_t)
#                    + b2 (PEPFAR_i x Post_t) [+ X_i x Post_t] + e_it
#
#   Y_it     Gross domestic savings (% of GDP), WDI NY.GDS.TOTL.ZS: GDP minus all final
#            consumption, so it excludes aid and remittances received from abroad
#   HIV_i    HIV prevalence ages 15-49 in 2003 (bl_hiv_prev), fixed before PEPFAR
#   PEPFAR_i 1 for the treatment group (All PEPFAR or COP), 0 for the Kates controls
#   Post_t   1 from 2004
#   HIV_i, PEPFAR_i and HIV_i x PEPFAR_i are absorbed by the country FE, Post_t by the
#   year FE, so only the three terms above are estimated.
#
#   b3 = how much more the savings rate rose after 2004 per percentage point of 2003 HIV
#        prevalence in PEPFAR countries than in control countries   <- main coefficient
#   b1 = the same HIV gradient among the controls
#   b1 + b3 = the HIV gradient among PEPFAR countries
#
# Adjusted adds Anthony's adjusted controls: the baseline (2003) covariates x Post and
# income group x Post, leaving out bl_hiv_prev, which is already in HIV_i x Post_t.
#
# Check. The controls have almost no HIV: in 2003, 60 of the 61 controls with HIV data
# were under 1% (median 0.1%). The exception, Equatorial Guinea (3.8%), has savings data
# only from 2005, so its country FE absorbs it and it does not help estimate b1. So b1
# rests on prevalence between 0 and 1%. The table reports:
#   Panel A  b3 from the full model
#   Panel B  b1 from the same regressions
#   Panel C  b3 from the model without HIV_i x Post_t, which sets the controls' HIV
#            gradient to zero, so b3 compares higher- and lower-HIV PEPFAR countries
#            with the control average
#
# Run from the repository root:  Rscript analysis/mechanisms/savings_panel_triple.R
# Needs data/intermediate/wdi_mechanisms_1990_2022.csv (data_pull/pull_wdi_mechanisms.R)
# Writes tables/crown_panel/mechanisms/savings_panel_triple.{csv,tex}

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/mechanisms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mech <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)

y_var      <- "savings_pct_gdp"
transforms <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs      <- c(unadj = "Unadjusted", adj = "Adjusted")

# Same transformations as fdi_panel.R (winsorizing uses percentiles of the full panel)
panel_m <- left_join(panel, mech, by = c("iso3c", "year"))
w <- winsor(panel_m[[y_var]])
panel_m[[paste0(y_var, "__raw")]]    <- panel_m[[y_var]]
panel_m[[paste0(y_var, "__win")]]    <- w
panel_m[[paste0(y_var, "__winihs")]] <- asinh(w)

# The three interaction terms. Countries without 2003 HIV data (3 controls) drop out,
# as they do in the HIV columns of fdi_panel.R.
make_triple_sample <- function(g) {
  make_sample(g, data = panel_m) |>
    filter(!is.na(bl_hiv_prev)) |>
    mutate(hiv_pepfar_post = bl_hiv_prev * treat * post,
           hiv_post        = bl_hiv_prev * post,
           pepfar_post     = treat * post)
}

# model: "full" (all three terms) or "no_hiv_post" (Panel C)
run_triple <- function(d, yvar, adjusted = FALSE, model = "full") {
  rhs <- c("hiv_pepfar_post", if (model == "full") "hiv_post", "pepfar_post")
  if (adjusted) {
    rhs <- c(rhs, paste0(setdiff(bl_cont, "bl_hiv_prev"), ":post"),
             "i(bl_wb_income, post, ref = 'L')")
  }
  feols(as.formula(paste(yvar, "~", paste(rhs, collapse = " + "), "| iso3c + year")),
        data = d, cluster = ~iso3c)
}

# One row per regression: b3, b1, b2 and b1 + b3, each with SE and p-value
tidy_triple <- function(m, d) {
  ct <- coeftable(m); b <- coef(m); V <- vcov(m)
  get <- function(term, what) if (term %in% rownames(ct)) ct[term, what] else NA_real_
  slope <- if ("hiv_post" %in% names(b)) {
    est <- b[["hiv_post"]] + b[["hiv_pepfar_post"]]
    se  <- sqrt(V["hiv_post", "hiv_post"] + V["hiv_pepfar_post", "hiv_pepfar_post"] +
                  2 * V["hiv_post", "hiv_pepfar_post"])
    # same t distribution fixest uses with country-clustered SEs (G - 1 df)
    c(est = est, se = se, p = 2 * pt(-abs(est / se), df = n_countries(m, d) - 1))
  } else c(est = NA_real_, se = NA_real_, p = NA_real_)
  tibble(
    b3 = get("hiv_pepfar_post", "Estimate"), se3 = get("hiv_pepfar_post", "Std. Error"),
    p3 = get("hiv_pepfar_post", "Pr(>|t|)"),
    b1 = get("hiv_post", "Estimate"), se1 = get("hiv_post", "Std. Error"),
    p1 = get("hiv_post", "Pr(>|t|)"),
    b2 = get("pepfar_post", "Estimate"), se2 = get("pepfar_post", "Std. Error"),
    p2 = get("pepfar_post", "Pr(>|t|)"),
    slope_pepfar = slope[["est"]], slope_pepfar_se = slope[["se"]],
    slope_pepfar_p = slope[["p"]],
    n_obs = nobs(m), n_countries = n_countries(m, d))
}

# =====================================================
# Estimate
# =====================================================
res <- list()
for (g in groups) {
  d <- make_triple_sample(g)
  for (sp in names(specs)) for (tr in names(transforms)) {
    yv  <- paste0(y_var, "__", tr)
    adj <- sp == "adj"
    fits <- list(full        = run_triple(d, yv, adj),
                 no_hiv_post = run_triple(d, yv, adj, model = "no_hiv_post"))
    for (mod in names(fits)) {
      res[[length(res) + 1]] <- bind_cols(
        tibble(outcome = y_var, group = g, model = mod, transform = tr, spec = sp),
        tidy_triple(fits[[mod]], d))
    }
  }
}
res <- bind_rows(res)
write_csv(res, file.path(out_dir, "savings_panel_triple.csv"))

# =====================================================
# LaTeX: one table. Columns: group x spec. Panels: A-C as described at the top.
# =====================================================
cols <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
pick <- function(mod, tr, j) filter(res, model == mod, transform == tr,
                                    group == cols$group[j], spec == cols$spec[j])
panel_rows <- function(mod, b = "b3", se = "se3", p = "p3", countries = TRUE) {
  out <- list()
  for (tr in names(transforms)) {
    r <- lapply(seq_len(nrow(cols)), \(j) pick(mod, tr, j))
    out[[length(out) + 1]] <- c(transforms[[tr]], sapply(r, \(x) tex_coef(x[[b]], x[[p]])))
    out[[length(out) + 1]] <- c("", sapply(r, \(x) tex_se(x[[se]])))
  }
  if (countries) {
    out[[length(out) + 1]] <- c("Countries", sapply(seq_len(nrow(cols)), \(j) {
      v <- filter(res, model == mod, group == cols$group[j], spec == cols$spec[j])$n_countries
      if (min(v) == max(v)) as.character(v[1]) else paste0(min(v), "--", max(v))
    }))
  }
  out
}
title_row <- function(txt) sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", txt)

rows <- c(
  list(title_row("Panel A: HIV $\\times$ PEPFAR $\\times$ Post")),
  panel_rows("full", countries = FALSE),
  list("\\midrule", title_row("Panel B: HIV $\\times$ Post (comparison countries), same regressions as A")),
  panel_rows("full", b = "b1", se = "se1", p = "p1"),
  list("\\midrule", title_row("Panel C: HIV $\\times$ PEPFAR $\\times$ Post, without HIV $\\times$ Post")),
  panel_rows("no_hiv_post"))

write_tex_table(
  file    = file.path(out_dir, "savings_panel_triple.tex"),
  caption = paste("Gross domestic savings (\\% of GDP): HIV prevalence $\\times$ PEPFAR",
                  "$\\times$ Post, two-way fixed effects"),
  label   = "tab:crown_mech_savings_triple",
  head    = list(
    c("", "\\multicolumn{2}{c}{All PEPFAR sample}", "\\multicolumn{2}{c}{COP sample}"),
    "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
    c("", rep(c("Unadj.", "Adj."), 2)),
    c("", sprintf("(%d)", 1:4))),
  rows    = rows,
  notes   = paste(
    "Two-way fixed effects regressions of gross domestic savings (\\% of GDP; GDP minus",
    "final consumption, World Bank WDI) with country and year fixed effects; Post = 1 from",
    "2004 for every country. HIV is HIV prevalence among ages 15--49 in 2003 (percentage",
    "points); PEPFAR is 1 for the treatment group of the column's sample (all PEPFAR",
    "recipients or COP countries) and 0 for the 64 comparison countries valid through 2022,",
    "of which the 61 with 2003 HIV data are used.",
    "Panel A: coefficient on HIV $\\times$ PEPFAR $\\times$ Post from regressions that also",
    "include HIV $\\times$ Post and PEPFAR $\\times$ Post. Panel B: coefficient on HIV",
    "$\\times$ Post from the same regressions, the HIV gradient among comparison countries.",
    "In 2003, 60 of the 61 comparison countries had HIV prevalence under 1\\% (median",
    "0.1\\%); the exception, Equatorial Guinea (3.8\\%), has savings data only from 2005.",
    "Panel C: Panel A without HIV $\\times$ Post, so the comparison countries' HIV gradient",
    "is set to zero.",
    "Adjusted columns add the baseline (2003) covariates of Crown et al.\\ (2023) other than",
    "HIV prevalence, each interacted with Post, and World Bank income group $\\times$ Post.",
    "Winsorized: values outside the 1st--99th percentiles of the full panel set to those",
    "percentiles. IHS: $\\sinh^{-1}(w)$ of the winsorized value, defined for the negative",
    "savings rates some countries record.",
    "Standard errors clustered by country in parentheses.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccc", resize = FALSE, size = "\\small")

cat("Savings triple-interaction table done. Output in", out_dir, "\n")

# Linear trend-break model for the mortality outcomes (Kates et al. replication)
# Author: Rishad Vaghaiwalla
# Follows Burgess and Pande (2005, AER), equation (3) and Table 1, with one break in 2004.
#
# Why: the mortality event studies show PEPFAR-vs-control gaps that drift in different
# directions before and after 2004, so a single before/after DiD can hide the change.
# This model fits a separate linear trend on each side of 2004 and tests the change.
#
# Model (one per outcome x treatment version x specification):
#   Y_it = a_i + g_t + γ1 (X_i x [t - 1990])
#                    + γ2 (X_i x [t - 2004] x Post_t)
#                    + γ4 (X_i x Post_t) + e_it
#   a_i = country FE, g_t = year FE, Post_t = 1 from 2004 on.
#   γ1 = pre-2004 trend of X countries relative to controls (per year)
#   γ2 = change in that trend from 2004 on   <- the coefficient we care about
#   γ4 = level jump in 2004, on top of the pre-2004 trend
#   F-test γ1 + γ2 = 0: is the post-2004 trend itself zero?
#   (γ4 keeps the name used in the assignment notes; there is no γ3 with one break.)
#
# Treatment versions X_i:
#   All PEPFAR  1 for any PEPFAR recipient, 0 for the Kates controls
#   COP         1 for COP countries, 0 for the Kates controls
#   HIV burden  HIV prevalence ages 15-49 in 2003 (bl_hiv_prev, % of population, so already
#               per-population), All PEPFAR + control sample. Fixed at its pre-PEPFAR value,
#               like Burgess and Pande's 1961 bank branches, so PEPFAR's own effect on
#               prevalence cannot feed back into X_i.
# Specifications:
#   Unadjusted  country and year FE only
#   Adjusted    the same controls as the adjusted event studies: the 11 baseline (2003)
#               covariates each with its own coefficient in every year, and income group x
#               year FE. In the HIV version bl_hiv_prev is left out of the controls, because
#               bl_hiv_prev x year dummies would absorb every bl_hiv_prev x trend term.
# Standard errors clustered by country.
#
# Run from the repository root:  Rscript analysis/parallel_trends/trend_break_mortality.R
# Writes to tables/kates_trend_breaks/{tex, figures, csv}/

source("analysis/econ_setup.R")

outcomes <- c(mort_allcause = "All-cause mortality",
              mort_child    = "Child mortality",
              mort_maternal = "Maternal mortality")
outcomes_tex <- c(mort_allcause = "All-cause",
                  mort_child    = "Child",
                  mort_maternal = "Maternal")

start_year <- 1990
break_year <- 2004   # Post_t = 1 from this year on (same as `post` in econ_setup.R)

out_dir <- "tables/kates_trend_breaks"
tex_dir <- file.path(out_dir, "tex")
fig_dir <- file.path(out_dir, "figures")
csv_dir <- file.path(out_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Treatment versions and model helpers
# =====================================================

versions <- list(
  "All PEPFAR" = list(data = samples[["All PEPFAR"]], xvar = "treat",
                      file = "all_pepfar", x_tex = "\\mathrm{PEPFAR}_i"),
  "COP"        = list(data = samples[["COP"]], xvar = "treat",
                      file = "cop", x_tex = "\\mathrm{COP}_i"),
  "HIV burden" = list(data = samples[["All PEPFAR"]], xvar = "bl_hiv_prev",
                      file = "hiv_burden", x_tex = "\\mathrm{HIV}_{i,2003}")
)
specs <- c("Unadjusted", "Adjusted")

# X_i and the three trend-break regressors
add_terms <- function(d, xvar) {
  d |>
    filter(!is.na(.data[[xvar]])) |>
    mutate(x         = .data[[xvar]],
           tr_base   = x * (year - start_year),           # γ1
           tr_change = x * (year - break_year) * post,    # γ2
           jump      = x * post)                          # γ4
}

# Fixed effects / controls part of the formula
fe_part <- function(spec, xvar) {
  if (spec == "Unadjusted") return("iso3c + year")
  ctrl <- if (xvar == "bl_hiv_prev") setdiff(bl_cont, "bl_hiv_prev") else bl_cont
  paste0("iso3c + bl_wb_income^year + year[", paste(ctrl, collapse = ", "), "]")
}

run_tb <- function(d, y, fe) {
  feols(as.formula(paste0(y, " ~ tr_base + tr_change + jump | ", fe)),
        data = d, cluster = ~iso3c)
}
# Event study with the same X_i, sample and controls, for the fit-check graphs
run_es_x <- function(d, y, fe) {
  feols(as.formula(paste0(y, " ~ i(year, x, ref = ", ref_year, ") | ", fe)),
        data = d, cluster = ~iso3c)
}

# =====================================================
# 2. Estimate
# =====================================================

res <- list()   # one row per trend-break model
es  <- list()   # event-study coefficients, for the graphs

for (v in names(versions)) for (s in specs) for (y in names(outcomes)) {
  vv <- versions[[v]]
  d  <- add_terms(vv$data, vv$xvar)
  fe <- fe_part(s, vv$xvar)

  m  <- run_tb(d, y, fe)
  b  <- coef(m); V <- vcov(m); ct <- coeftable(m)
  G  <- n_countries(m, d)

  # γ1 + γ2 (the post-2004 trend) and its F-test, df = (1, G - 1) as fixest uses for
  # country-clustered SEs
  w   <- c(tr_base = 1, tr_change = 1, jump = 0)[names(b)]
  sum12    <- sum(w * b)
  sum12_se <- sqrt(drop(t(w) %*% V %*% w))
  F12      <- (sum12 / sum12_se)^2

  used <- d[obs(m), ]
  res[[length(res) + 1]] <- tibble(
    version = v, spec = s, outcome = y,
    g1 = b[["tr_base"]],   g1_se = ct["tr_base", "Std. Error"],   g1_p = ct["tr_base", "Pr(>|t|)"],
    g2 = b[["tr_change"]], g2_se = ct["tr_change", "Std. Error"], g2_p = ct["tr_change", "Pr(>|t|)"],
    g4 = b[["jump"]],      g4_se = ct["jump", "Std. Error"],      g4_p = ct["jump", "Pr(>|t|)"],
    post_trend = sum12, post_trend_se = sum12_se,
    F_g1g2 = F12, F_g1g2_p = pf(F12, 1, G - 1, lower.tail = FALSE),
    adj_r2 = unname(r2(m, "ar2")),
    y_mean_2003 = mean(used[[y]][used$year == ref_year], na.rm = TRUE),
    n_obs = nobs(m), n_countries = G)

  mes <- run_es_x(d, y, fe)
  ect <- coeftable(mes)
  keep <- grepl("^year::\\d+:x$", rownames(ect))
  es[[length(es) + 1]] <- tibble(
    version = v, spec = s, outcome = y,
    year = as.integer(sub("year::(\\d+):x", "\\1", rownames(ect)[keep])),
    est = ect[keep, "Estimate"], se = ect[keep, "Std. Error"]) |>
    bind_rows(tibble(version = v, spec = s, outcome = y, year = ref_year, est = 0, se = 0))
}
res <- bind_rows(res)
es  <- bind_rows(es) |> arrange(version, spec, outcome, year) |>
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se)
write_csv(res, file.path(csv_dir, "trend_break_results.csv"))
write_csv(es,  file.path(csv_dir, "trend_break_event_study_coefficients.csv"))

# =====================================================
# 3. LaTeX tables (Burgess and Pande Table 1 layout), one per treatment version
#    Columns (1)-(3) unadjusted, (4)-(6) adjusted; each block all-cause, child, maternal
# =====================================================

ys   <- names(outcomes)
cols <- expand.grid(outcome = ys, spec = specs, stringsAsFactors = FALSE)

cell <- function(v, s, y, what) {
  r <- filter(res, version == v, spec == s, outcome == y)
  switch(what,
    g1 = tex_coef(r$g1, r$g1_p), g1_se = tex_se(r$g1_se),
    g2 = tex_coef(r$g2, r$g2_p), g2_se = tex_se(r$g2_se),
    g4 = tex_coef(r$g4, r$g4_p), g4_se = tex_se(r$g4_se),
    F  = sprintf("%.2f", r$F_g1g2),
    Fp = sprintf("[%.3f]", r$F_g1g2_p),
    ymean = sprintf("%.1f", r$y_mean_2003),
    r2 = sprintf("%.3f", r$adj_r2),
    n  = tex_n(r$n_obs), G = tex_n(r$n_countries))
}
row_of <- function(label, v, what) {
  c(label, mapply(function(s, y) cell(v, s, y, what), cols$spec, cols$outcome))
}

x_note <- c(
  "All PEPFAR" = "$\\mathrm{PEPFAR}_i$ is 1 for every PEPFAR recipient and 0 for the controls.",
  "COP"        = "$\\mathrm{COP}_i$ is 1 for COP countries and 0 for the controls.",
  "HIV burden" = paste(
    "$\\mathrm{HIV}_{i,2003}$ is HIV prevalence among ages 15--49 in 2003, in percent",
    "(WDI, gap-filled with GBD), so coefficients are per percentage point of 2003 prevalence;",
    "sample is all PEPFAR recipients and controls with 2003 prevalence.",
    "Baseline HIV prevalence is not among the adjusted controls, since its year-by-year",
    "coefficients would absorb the trend terms."))

for (v in names(versions)) {
  xt <- versions[[v]]$x_tex
  rows <- list(
    row_of(sprintf("$%s \\times (t-1990)$ \\quad [$\\gamma_1$]", xt), v, "g1"),
    row_of("", v, "g1_se"),
    row_of(sprintf("$%s \\times (t-2004) \\times \\mathrm{Post}_t$ \\quad [$\\gamma_2$]", xt), v, "g2"),
    row_of("", v, "g2_se"),
    row_of(sprintf("$%s \\times \\mathrm{Post}_t$ \\quad [$\\gamma_4$]", xt), v, "g4"),
    row_of("", v, "g4_se"),
    "\\midrule",
    c("Country and year FE", rep("Yes", 6)),
    c("Baseline covariates $\\times$ year", rep("No", 3), rep("Yes", 3)),
    row_of("F-test: $\\gamma_1 + \\gamma_2 = 0$", v, "F"),
    row_of("", v, "Fp"),
    row_of("Outcome mean, 2003", v, "ymean"),
    row_of("Adjusted $R^2$", v, "r2"),
    row_of("Observations", v, "n"),
    row_of("Countries", v, "G"))

  write_tex_table(
    file    = file.path(tex_dir, paste0("trend_break_", versions[[v]]$file, ".tex")),
    caption = paste0("Linear trend breaks in mortality around 2004: ",
                     if (v == "HIV burden") "baseline HIV prevalence" else paste(v, "vs.\\ controls")),
    label   = paste0("tab:trend_break_", versions[[v]]$file),
    head    = list(c("", "\\multicolumn{3}{c}{Unadjusted}", "\\multicolumn{3}{c}{Adjusted}"),
                   "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
                   c("", sprintf("(%d)", 1:6)),
                   c("", rep(outcomes_tex[ys], 2))),
    rows    = rows,
    notes   = paste(
      "Estimates of $Y_{it} = \\alpha_i + \\gamma_t + \\gamma_1 (X_i \\times [t-1990])",
      "+ \\gamma_2 (X_i \\times [t-2004] \\times \\mathrm{Post}_t) + \\gamma_4 (X_i \\times \\mathrm{Post}_t)",
      "+ \\varepsilon_{it}$, following Burgess and Pande (2005), with $\\mathrm{Post}_t = 1$ from 2004.",
      "$\\gamma_1$ is the pre-2004 trend relative to controls (per year), $\\gamma_2$ the change in that",
      "trend from 2004, and $\\gamma_4$ the level shift in 2004.",
      "The F-test (1, countries $-$ 1 df) is that the post-2004 trend $\\gamma_1 + \\gamma_2$ is zero;",
      "$p$-value in brackets.", x_note[[v]],
      if (v == "HIV burden") "Adjusted columns add the other 10 baseline (2003) covariates" else
        "Adjusted columns add the 11 baseline (2003) covariates",
      "of Crown et al.\\ (2023), each with its own coefficient in every year, and World Bank",
      "income group $\\times$ year fixed effects.",
      "Outcomes: crude death rate per 1,000 people, under-5 deaths per 1,000 live births, and",
      "maternal deaths per 100,000 live births (WDI).",
      "Sample 1990--2022 (maternal mortality from 2000); controls are the 64 comparison countries",
      "valid through 2022 (\\texttt{grp\\_control\\_kates}). Standard errors clustered by country",
      "in parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
    align   = "lcccccc")
}

# =====================================================
# 4. Fit check: event-study coefficients with the fitted trend-break line
#    Fitted gap vs. 2003 = γ1 (k - 2003) + 1[k >= 2004] (γ4 + γ2 (k - 2004))
# =====================================================

fit_line <- res |>
  select(version, spec, outcome, g1, g2, g4) |>
  inner_join(distinct(es, version, spec, outcome, year),
             by = c("version", "spec", "outcome"), relationship = "one-to-many") |>
  mutate(fit = g1 * (year - ref_year) +
               (year >= break_year) * (g4 + g2 * (year - break_year)),
         side = ifelse(year >= break_year, "post", "pre"))

ver_lab <- c("All PEPFAR" = "All PEPFAR", "COP" = "COP",
             "HIV burden" = "HIV prevalence 2003 (per pp)")

save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, bg = "white")
}

fit_plot <- function(s) {
  lab <- function(d) d |> mutate(
    panel = factor(paste0(outcomes[outcome], "\n", ver_lab[version]),
                   levels = as.vector(outer(ver_lab, outcomes, function(a, b) paste0(b, "\n", a)))))
  e <- lab(filter(es, spec == s))
  f <- lab(filter(fit_line, spec == s))
  ggplot(e, aes(year, est)) +
    geom_hline(yintercept = 0, color = ink2, linewidth = 0.4) +
    geom_vline(xintercept = break_year - 0.5, color = ink2, linetype = "dashed",
               linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, color = pal[1],
                  alpha = 0.6) +
    geom_point(size = 1.2, color = pal[1]) +
    geom_line(data = f, aes(year, fit, group = side), color = pal[2], linewidth = 0.9) +
    facet_wrap(~panel, scales = "free_y", ncol = 3) +
    labs(title = paste0(s, " event studies with fitted trend breaks"),
         subtitle = "Dots: event-study gap vs. 2003 (95% CI). Line: fitted trend-break model.",
         caption = paste0("Country and year FE",
                          if (s == "Adjusted") " + baseline covariates x year, income group x year FE",
                          " | SEs clustered by country | HIV panels: per percentage point of 2003 prevalence"),
         x = NULL, y = "Gap vs. 2003") +
    theme_es +
    theme(plot.caption = element_text(color = ink2, hjust = 0, size = 7),
          strip.text = element_text(size = 8))
}
for (s in specs) save_both(fit_plot(s), paste0("trend_break_fit_", tolower(s)), 9, 8.5)

# =====================================================
# 5. LaTeX figure blocks
# =====================================================
fig_env <- function(file, caption, label) {
  c("\\begin{figure}[htbp]", "\\centering",
    sprintf("\\includegraphics[width=\\textwidth]{figures/%s.pdf}", file),
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label), "\\end{figure}")
}
writeLines(c(
  fig_env("trend_break_fit_unadjusted",
          "Unadjusted mortality event studies with fitted linear trend breaks",
          "fig:trend_break_fit_unadjusted"),
  fig_env("trend_break_fit_adjusted",
          "Adjusted mortality event studies with fitted linear trend breaks",
          "fig:trend_break_fit_adjusted")),
  file.path(tex_dir, "figures_trend_break.tex"))

cat("Trend-break results done. Output in", out_dir, "\n")

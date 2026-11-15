# Maternal mortality trend breaks: 2004 vs. 2008, plus a break-year sweep
# Author: Rishad Vaghaiwalla
#
# Same model as trend_break_mortality.R (Burgess and Pande 2005, eq. 3), maternal mortality
# only, with the break year B as a parameter:
#   Y_it = a_i + g_t + γ1 (X_i x [t - 1990])
#                    + γ2 (X_i x [t - B] x Post_t(B))
#                    + γ4 (X_i x Post_t(B)) + e_it,     Post_t(B) = 1 from B on
# Treatment versions (All PEPFAR, COP, HIV burden) and specifications (Unadjusted,
# Adjusted) are exactly as in trend_break_mortality.R, so the B = 2004 columns reproduce
# tables/kates_trend_breaks/csv/trend_break_results.csv.
#
# Why 2008: proposed in the team meeting for maternal mortality. Maternal mortality data
# start in 2000, so with B = 2004 the pre-break slope γ1 rests on only 2000-2003.
# The sweep refits every version x spec for each break year from 2002 to 2016 and records
# γ2 and the residual sum of squares (RSS). The break year with the lowest RSS is the
# least-squares estimate of the break date; its p-values are not valid as-is, because
# the data chose the year.
#
# Run from the repository root:  Rscript analysis/parallel_trends/trend_break_maternal_2008.R
# Writes to tables/kates_trend_breaks/maternal_2008/{tex, figures, csv}/

source("analysis/econ_setup.R")

y          <- "mort_maternal"
start_year <- 1990
breaks_cmp <- c(2004, 2008)   # the two break years compared side by side
sweep_yrs  <- 2002:2016       # candidate break years for the sweep

out_dir <- "tables/kates_trend_breaks/maternal_2008"
tex_dir <- file.path(out_dir, "tex")
fig_dir <- file.path(out_dir, "figures")
csv_dir <- file.path(out_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Model helpers (same as trend_break_mortality.R, with B as an argument)
# =====================================================

versions <- list(
  "All PEPFAR" = list(data = samples[["All PEPFAR"]], xvar = "treat"),
  "COP"        = list(data = samples[["COP"]],        xvar = "treat"),
  "HIV burden" = list(data = samples[["All PEPFAR"]], xvar = "bl_hiv_prev")
)
specs <- c("Unadjusted", "Adjusted")

add_terms <- function(d, xvar, B) {
  d |>
    filter(!is.na(.data[[xvar]])) |>
    mutate(x         = .data[[xvar]],
           post_b    = as.integer(year >= B),
           tr_base   = x * (year - start_year),     # γ1
           tr_change = x * (year - B) * post_b,     # γ2
           jump      = x * post_b)                  # γ4
}

fe_part <- function(spec, xvar) {
  if (spec == "Unadjusted") return("iso3c + year")
  ctrl <- if (xvar == "bl_hiv_prev") setdiff(bl_cont, "bl_hiv_prev") else bl_cont
  paste0("iso3c + bl_wb_income^year + year[", paste(ctrl, collapse = ", "), "]")
}

run_tb <- function(d, fe) {
  feols(as.formula(paste0(y, " ~ tr_base + tr_change + jump | ", fe)),
        data = d, cluster = ~iso3c)
}

# One row of results for version v, spec s, break year B
fit_row <- function(v, s, B) {
  vv <- versions[[v]]
  d  <- add_terms(vv$data, vv$xvar, B)
  m  <- run_tb(d, fe_part(s, vv$xvar))
  b  <- coef(m); V <- vcov(m); ct <- coeftable(m)
  G  <- n_countries(m, d)
  w  <- c(tr_base = 1, tr_change = 1, jump = 0)[names(b)]
  s12 <- sum(w * b); s12_se <- sqrt(drop(t(w) %*% V %*% w))
  F12 <- (s12 / s12_se)^2
  tibble(version = v, spec = s, break_year = B,
         g1 = b[["tr_base"]],   g1_se = ct["tr_base", "Std. Error"],   g1_p = ct["tr_base", "Pr(>|t|)"],
         g2 = b[["tr_change"]], g2_se = ct["tr_change", "Std. Error"], g2_p = ct["tr_change", "Pr(>|t|)"],
         g4 = b[["jump"]],      g4_se = ct["jump", "Std. Error"],      g4_p = ct["jump", "Pr(>|t|)"],
         post_trend = s12, post_trend_se = s12_se,
         F_g1g2 = F12, F_g1g2_p = pf(F12, 1, G - 1, lower.tail = FALSE),
         rss = sum(resid(m)^2), adj_r2 = unname(r2(m, "ar2")),
         n_obs = nobs(m), n_countries = G)
}

# =====================================================
# 2. 2004 vs. 2008
# =====================================================

cmp <- bind_rows(lapply(breaks_cmp, function(B)
  bind_rows(lapply(names(versions), function(v)
    bind_rows(lapply(specs, function(s) fit_row(v, s, B)))))))
write_csv(cmp, file.path(csv_dir, "maternal_break_2004_vs_2008.csv"))

# LaTeX: columns = version x spec, Panel A = 2004, Panel B = 2008
cols <- expand.grid(spec = specs, version = names(versions), stringsAsFactors = FALSE)
cell <- function(B, v, s, what) {
  r <- filter(cmp, break_year == B, version == v, spec == s)
  switch(what,
         g1 = tex_coef(r$g1, r$g1_p), g1_se = tex_se(r$g1_se),
         g2 = tex_coef(r$g2, r$g2_p), g2_se = tex_se(r$g2_se),
         g4 = tex_coef(r$g4, r$g4_p), g4_se = tex_se(r$g4_se),
         F  = sprintf("%.2f", r$F_g1g2), Fp = sprintf("[%.3f]", r$F_g1g2_p),
         n  = tex_n(r$n_obs), G = tex_n(r$n_countries))
}
row_of <- function(label, B, what) c(label, mapply(function(v, s) cell(B, v, s, what),
                                                   cols$version, cols$spec))
panel_rows <- function(B) list(
  sprintf("\\multicolumn{7}{l}{\\textit{Break in %d}} \\\\", B),
  row_of("$X_i \\times (t-1990)$ \\quad [$\\gamma_1$]", B, "g1"),     row_of("", B, "g1_se"),
  row_of("$X_i \\times (t-B) \\times \\mathrm{Post}_t$ \\quad [$\\gamma_2$]", B, "g2"), row_of("", B, "g2_se"),
  row_of("$X_i \\times \\mathrm{Post}_t$ \\quad [$\\gamma_4$]", B, "g4"), row_of("", B, "g4_se"),
  row_of("F-test: $\\gamma_1 + \\gamma_2 = 0$", B, "F"), row_of("", B, "Fp"),
  row_of("Observations", B, "n"), row_of("Countries", B, "G"))

write_tex_table(
  file    = file.path(tex_dir, "trend_break_maternal_2004_vs_2008.tex"),
  caption = "Linear trend breaks in maternal mortality: break in 2004 vs.\\ 2008",
  label   = "tab:trend_break_maternal_2008",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR}", "\\multicolumn{2}{c}{COP}",
                   "\\multicolumn{2}{c}{HIV prevalence 2003}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
                 c("", rep(c("Unadj.", "Adj."), 3)),
                 c("", sprintf("(%d)", 1:6))),
  rows    = c(panel_rows(2004), list("\\midrule"), panel_rows(2008)),
  notes   = paste(
    "Estimates of $Y_{it} = \\alpha_i + \\gamma_t + \\gamma_1 (X_i \\times [t-1990])",
    "+ \\gamma_2 (X_i \\times [t-B] \\times \\mathrm{Post}_t) + \\gamma_4 (X_i \\times \\mathrm{Post}_t)",
    "+ \\varepsilon_{it}$ with $\\mathrm{Post}_t = 1$ from the break year $B$.",
    "Outcome: maternal deaths per 100,000 live births (WDI), 2000--2022.",
    "$X_i$ is the PEPFAR or COP indicator, or HIV prevalence among ages 15--49 in 2003 (per",
    "percentage point). Adjusted columns add the baseline (2003) covariates of Crown et al.\\ (2023),",
    "each with its own coefficient in every year, and World Bank income group $\\times$ year fixed",
    "effects (baseline HIV prevalence is left out of the controls in columns 5--6).",
    "The F-test (1, countries $-$ 1 df) is that the post-break trend $\\gamma_1 + \\gamma_2$ is zero;",
    "$p$-value in brackets. Standard errors clustered by country in parentheses.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccccc")

# =====================================================
# 3. Break-year sweep
# =====================================================

sweep <- bind_rows(lapply(sweep_yrs, function(B)
  bind_rows(lapply(names(versions), function(v)
    bind_rows(lapply(specs, function(s) fit_row(v, s, B))))))) |>
  group_by(version, spec) |>
  mutate(rss_rel  = rss / min(rss),                  # 1 = best-fitting break year
         best_rss = break_year == break_year[which.min(rss)]) |>
  ungroup()
write_csv(sweep, file.path(csv_dir, "maternal_break_sweep.csv"))

# Same sample for every break year? (otherwise RSS is not comparable)
stopifnot(sweep |> group_by(version, spec) |> summarise(k = n_distinct(n_obs), .groups = "drop") |>
            pull(k) |> max() == 1)

best <- sweep |> filter(best_rss) |> select(version, spec, best_break = break_year)
write_csv(best, file.path(csv_dir, "maternal_best_break_year.csv"))

save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
sw <- sweep |> mutate(panel = factor(paste(version, spec, sep = " | "),
                                     levels = paste(rep(names(versions), each = 2), specs, sep = " | ")),
                      lo = g2 - 1.96 * g2_se, hi = g2 + 1.96 * g2_se)

p_g2 <- ggplot(sw, aes(break_year, g2)) +
  geom_hline(yintercept = 0, color = ink2, linewidth = 0.4) +
  geom_vline(xintercept = breaks_cmp, color = ink2, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, color = pal[1], alpha = 0.6) +
  geom_point(color = pal[1], size = 1.4) +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(title = "Maternal mortality: slope change (gamma 2) by break year",
       subtitle = "95% CI, SEs clustered by country. Dashed lines: 2004 and 2008.",
       x = "Break year B", y = "Slope change, gamma 2") +
  theme_es
save_both(p_g2, "maternal_sweep_g2", 8, 8)

p_rss <- ggplot(sw, aes(break_year, rss_rel)) +
  geom_vline(xintercept = breaks_cmp, color = ink2, linetype = "dashed", linewidth = 0.3) +
  geom_line(color = pal[2]) + geom_point(aes(shape = best_rss), color = pal[2], size = 1.6) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), labels = c("", "Lowest RSS")) +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(title = "Maternal mortality: model fit by break year",
       subtitle = "Residual sum of squares relative to the best break year (1 = best fit)",
       x = "Break year B", y = "RSS / min RSS") +
  theme_es + theme(legend.position = "none")
save_both(p_rss, "maternal_sweep_rss", 8, 8)

# =====================================================
# 4. Fit check: event study with both fitted trend-break lines
#    Fitted gap vs. 2003 = γ1 (k - 2003) + 1[k >= B] (γ4 + γ2 (k - B))
# =====================================================

es <- bind_rows(lapply(names(versions), function(v) bind_rows(lapply(specs, function(s) {
  vv  <- versions[[v]]
  d   <- add_terms(vv$data, vv$xvar, 2004)
  m   <- feols(as.formula(paste0(y, " ~ i(year, x, ref = ", ref_year, ") | ", fe_part(s, vv$xvar))),
               data = d, cluster = ~iso3c)
  ct  <- coeftable(m)
  keep <- grepl("^year::\\d+:x$", rownames(ct))
  tibble(version = v, spec = s,
         year = as.integer(sub("year::(\\d+):x", "\\1", rownames(ct)[keep])),
         est = ct[keep, "Estimate"], se = ct[keep, "Std. Error"]) |>
    bind_rows(tibble(version = v, spec = s, year = ref_year, est = 0, se = 0))
})))) |> mutate(lo = est - 1.96 * se, hi = est + 1.96 * se)

fits <- cmp |> select(version, spec, break_year, g1, g2, g4) |>
  inner_join(distinct(es, version, spec, year), by = c("version", "spec"),
             relationship = "many-to-many") |>
  mutate(fit  = g1 * (year - ref_year) + (year >= break_year) * (g4 + g2 * (year - break_year)),
         side = ifelse(year >= break_year, "post", "pre"),
         Break = factor(break_year))
lab <- function(d) d |> mutate(panel = factor(paste(version, spec, sep = " | "), levels = levels(sw$panel)))

p_fit <- ggplot(lab(es), aes(year, est)) +
  geom_hline(yintercept = 0, color = ink2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, color = ink2, alpha = 0.5) +
  geom_point(size = 1.1, color = ink2) +
  geom_line(data = lab(fits), aes(year, fit, color = Break, group = interaction(Break, side)),
            linewidth = 0.9) +
  scale_color_manual(values = pal[1:2], name = "Break year") +
  facet_wrap(~panel, scales = "free_y", ncol = 2) +
  labs(title = "Maternal mortality event studies with fitted trend breaks",
       subtitle = "Dots: event-study gap vs. 2003 (95% CI). Lines: fitted model with a break in 2004 or 2008.",
       x = NULL, y = "Gap vs. 2003 (deaths per 100,000 live births)") +
  theme_es + theme(legend.title = element_text(color = ink))
save_both(p_fit, "maternal_fit_2004_vs_2008", 8, 8.5)

# =====================================================
# 5. Console summary
# =====================================================
cat("\n2004 vs. 2008 (maternal mortality):\n")
print(cmp |> select(break_year, version, spec, g1, g2, g4, F_g1g2_p, n_obs) |>
        mutate(across(c(g1, g2, g4, F_g1g2_p), ~ round(.x, 3))), n = Inf)
cat("\nLowest-RSS break year by version x spec:\n")
print(best, n = Inf)
cat("\nDone. Output in", out_dir, "\n")
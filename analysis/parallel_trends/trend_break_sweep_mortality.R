# Step 4: break-year sweep for child, all-cause and maternal mortality
# Author: Rishad Vaghaiwalla
#
# trend_break_maternal_2008.R, extended from maternal mortality to all three outcomes.
# Same model and helpers (Burgess and Pande 2005, eq. 3), with the break year B as a parameter:
#   Y_it = a_i + g_t + γ1 (X_i x [t - 1990]) + γ2 (X_i x [t - B] x Post_t(B))
#                    + γ4 (X_i x Post_t(B)) + e_it
# For every outcome x treatment version x specification it refits the model for each B from
# 2002 to 2016 and records γ2, γ4, the residual sum of squares (RSS) and the joint Wald test
# of γ2 = γ4 = 0 (no break at all).
#
# How to pick the break year:
#   - The B with the lowest RSS is the least-squares estimate of the break date (Bai 1997).
#   - The largest Wald statistic over B is the Quandt-Andrews sup-Wald statistic. Its
#     p-value is NOT the usual one, because the data chose B; compare it with the
#     Andrews (1993) critical values for 2 restrictions and the trimming used here.
#   - The break used downstream is restricted to B >= 2004: a break before PEPFAR existed
#     cannot be credited to PEPFAR. The unrestricted best year is reported alongside.
#     The fit surface is often flat (RSS within a fraction of a percent across several
#     years), so report how much worse 2004 and 2008 fit, not just the winner.
#   - The year should also match a program milestone (gameplan Step 4 table). The maternal
#     rows reproduce tables/kates_trend_breaks/maternal_2008/csv/maternal_break_sweep.csv.
#
# Writes best_break_year.csv, which life_years_mortality.R and cost_benefit_income_groups.R
# read to set the child and all-cause break years (HIV burden version, adjusted model).
#
# Run from the repository root:  Rscript analysis/parallel_trends/trend_break_sweep_mortality.R
# Writes to tables/kates_trend_breaks/sweep/{tex, figures, csv}/

source("analysis/econ_setup.R")

outcomes <- c(mort_child    = "Child mortality (under 5)",
              mort_allcause = "All-cause mortality",
              mort_maternal = "Maternal mortality")
start_year <- 1990
sweep_yrs  <- 2002:2016       # candidate break years (same window as the maternal sweep)
min_break  <- 2004            # the break used downstream must be on or after PEPFAR's first year
specs      <- c("Unadjusted", "Adjusted")

out_dir <- "tables/kates_trend_breaks/sweep"
tex_dir <- file.path(out_dir, "tex"); fig_dir <- file.path(out_dir, "figures")
csv_dir <- file.path(out_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Model helpers (from trend_break_maternal_2008.R, with the outcome as an argument)
# =====================================================
versions <- list(
  "All PEPFAR" = list(data = samples[["All PEPFAR"]], xvar = "treat"),
  "COP"        = list(data = samples[["COP"]],        xvar = "treat"),
  "HIV burden" = list(data = samples[["All PEPFAR"]], xvar = "bl_hiv_prev")
)

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

fit_row <- function(y, v, s, B) {
  vv <- versions[[v]]
  d  <- add_terms(vv$data, vv$xvar, B)
  m  <- feols(as.formula(paste0(y, " ~ tr_base + tr_change + jump | ", fe_part(s, vv$xvar))),
              data = d, cluster = ~iso3c)
  b  <- coef(m); ct <- coeftable(m)
  wd <- wald(m, keep = "^tr_change$|^jump$", print = FALSE)   # H0: γ2 = γ4 = 0
  tibble(outcome = y, version = v, spec = s, break_year = B,
         g1 = b[["tr_base"]],   g1_se = ct["tr_base", "Std. Error"],
         g2 = b[["tr_change"]], g2_se = ct["tr_change", "Std. Error"], g2_p = ct["tr_change", "Pr(>|t|)"],
         g4 = b[["jump"]],      g4_se = ct["jump", "Std. Error"],      g4_p = ct["jump", "Pr(>|t|)"],
         wald_F = unname(wd$stat), wald_p_naive = unname(wd$p),
         rss = sum(resid(m)^2), n_obs = nobs(m), n_countries = n_countries(m, d))
}

# =====================================================
# 2. Sweep
# =====================================================
sweep <- bind_rows(lapply(names(outcomes), function(y)
  bind_rows(lapply(sweep_yrs, function(B)
    bind_rows(lapply(names(versions), function(v)
      bind_rows(lapply(specs, function(s) fit_row(y, v, s, B))))))))) |>
  group_by(outcome, version, spec) |>
  mutate(rss_rel  = rss / min(rss),                         # 1 = best-fitting break year
         pct_above_best = 100 * (rss_rel - 1),              # % worse fit than the best year
         best_rss = break_year == break_year[which.min(rss)],
         sup_wald = wald_F == max(wald_F)) |>
  ungroup()
write_csv(sweep, file.path(csv_dir, "break_sweep.csv"))

# Same sample for every break year? (otherwise RSS is not comparable across B)
stopifnot(sweep |> group_by(outcome, version, spec) |>
            summarise(k = n_distinct(n_obs), .groups = "drop") |> pull(k) |> max() == 1)

best <- sweep |>
  group_by(outcome, version, spec) |>
  summarise(best_break_any = break_year[which.min(rss)],               # unrestricted
            best_break     = break_year[break_year >= min_break][which.min(rss[break_year >= min_break])],
            sup_wald_year = break_year[which.max(wald_F)],
            sup_wald_F    = max(wald_F),
            g2_at_best    = g2[which.min(rss)], g2_se_at_best = g2_se[which.min(rss)],
            pct_above_best_2004 = pct_above_best[break_year == 2004],
            pct_above_best_2008 = pct_above_best[break_year == 2008],
            .groups = "drop")
write_csv(best, file.path(csv_dir, "best_break_year.csv"))

# =====================================================
# 3. Figures: RSS and γ2 by break year, one file per outcome (6 panels: version x spec)
# =====================================================
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
panel_levels <- paste(rep(names(versions), each = 2), specs, sep = " | ")
sw <- sweep |> mutate(panel = factor(paste(version, spec, sep = " | "), levels = panel_levels),
                      lo = g2 - 1.96 * g2_se, hi = g2 + 1.96 * g2_se)
ofile <- c(mort_child = "child", mort_allcause = "allcause", mort_maternal = "maternal")

for (y in names(outcomes)) {
  d <- filter(sw, outcome == y)
  p_rss <- ggplot(d, aes(break_year, pct_above_best)) +
    geom_vline(xintercept = c(2004, 2008), color = ink2, linetype = "dashed", linewidth = 0.3) +
    geom_line(color = pal[2]) + geom_point(aes(shape = best_rss), color = pal[2], size = 1.6) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8)) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    labs(title = paste0(outcomes[[y]], ": model fit by break year"),
         subtitle = "RSS, % above the best break year (0 = best fit, star). Dashed: 2004 and 2008.",
         x = "Break year B", y = "RSS, % above best") +
    theme_es + theme(legend.position = "none")
  save_both(p_rss, paste0(ofile[[y]], "_sweep_rss"), 8, 8)

  p_g2 <- ggplot(d, aes(break_year, g2)) +
    geom_hline(yintercept = 0, color = ink2, linewidth = 0.4) +
    geom_vline(xintercept = c(2004, 2008), color = ink2, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, color = pal[1], alpha = 0.6) +
    geom_point(color = pal[1], size = 1.4) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    labs(title = paste0(outcomes[[y]], ": slope change (gamma 2) by break year"),
         subtitle = "95% CI, SEs clustered by country. Dashed lines: 2004 and 2008.",
         x = "Break year B", y = "Slope change, gamma 2") +
    theme_es
  save_both(p_g2, paste0(ofile[[y]], "_sweep_g2"), 8, 8)
}

# =====================================================
# 4. LaTeX: best break year, sup-Wald year and statistic, for every outcome x version x spec
# =====================================================
cols <- expand.grid(spec = specs, version = names(versions), stringsAsFactors = FALSE)
cell <- function(y, v, s, what) {
  r <- filter(best, outcome == y, version == v, spec == s)
  switch(what, best = as.character(r$best_break), best_any = as.character(r$best_break_any),
         sup = sprintf("%d (%.1f)", r$sup_wald_year, r$sup_wald_F),
         r04 = sprintf("%.3f", r$pct_above_best_2004), r08 = sprintf("%.3f", r$pct_above_best_2008))
}
row_of <- function(label, y, what) c(label, mapply(function(v, s) cell(y, v, s, what), cols$version, cols$spec))
rows <- unlist(lapply(names(outcomes), function(y) list(
  sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", outcomes[[y]]),
  row_of("Lowest-RSS break year, 2004 on", y, "best"),
  row_of("Lowest-RSS break year, any", y, "best_any"),
  row_of("sup-Wald year (F)", y, "sup"),
  row_of("RSS \\% above best: 2004", y, "r04"),
  row_of("RSS \\% above best: 2008", y, "r08"),
  "\\midrule")), recursive = FALSE)
rows <- rows[-length(rows)]

write_tex_table(
  file    = file.path(tex_dir, "break_year_sweep.tex"),
  caption = "Best-fitting trend-break year by mortality outcome",
  label   = "tab:break_year_sweep",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR}", "\\multicolumn{2}{c}{COP}",
                   "\\multicolumn{2}{c}{HIV prevalence 2003}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
                 c("", rep(c("Unadj.", "Adj."), 3))),
  rows    = rows,
  notes   = paste(
    "Each cell refits $Y_{it} = \\alpha_i + \\gamma_t + \\gamma_1 (X_i \\times [t-1990])",
    "+ \\gamma_2 (X_i \\times [t-B] \\times \\mathrm{Post}_t) + \\gamma_4 (X_i \\times \\mathrm{Post}_t)",
    "+ \\varepsilon_{it}$ for every break year $B$ from", min(sweep_yrs), "to", paste0(max(sweep_yrs), "."),
    "Lowest-RSS year: least-squares estimate of the break date, restricted to", min_break,
    "or later (PEPFAR's first year) and unrestricted. sup-Wald: year with the largest",
    "joint Wald statistic for $\\gamma_2 = \\gamma_4 = 0$ (statistic in parentheses); compare with",
    "Andrews (1993) critical values, not the usual F distribution, since the break year is chosen",
    "by the data. RSS rows: percent by which the fit at 2004 and 2008 is worse than the best year.",
    "Samples and controls as in the trend-break tables; standard errors clustered by country."),
  align   = "lcccccc")

# =====================================================
# 5. Console summary
# =====================================================
cat("\nBest break year (lowest RSS) and sup-Wald year, by outcome x version x spec:\n")
print(best |> mutate(across(c(sup_wald_F, g2_at_best, pct_above_best_2004, pct_above_best_2008), ~ round(.x, 3))),
      n = Inf, width = Inf)
cat("\nUsed by the life-years scripts (HIV burden, Adjusted):\n")
print(filter(best, version == "HIV burden", spec == "Adjusted") |>
        select(outcome, best_break, best_break_any, pct_above_best_2004, pct_above_best_2008))
cat("(maternal uses 2008 from the team meeting, set in life_years_setup.R)\n")
cat("\nDone. Output in", out_dir, "\n")

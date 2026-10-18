# Event studies for the Crown et al. (2023) economic outcomes, 1990-2022
# Author: Anthony Babu
#
# Why: the 2x2 DiD (did_crown.R) gives one number, PEPFAR x Post, which lumps all
# pre years together and all post years together. An event study estimates a
# separate PEPFAR-vs-control gap for every year, relative to 2003, so we can see
# (1) whether the groups were already diverging before 2004 (parallel trends) and
# (2) when, after 2004, any effect shows up.
#
# Main model (one per outcome x treatment group x spec):
#   Y_it = a_i + g_t + sum_{k != 2003} b_k * (PEPFAR_i x 1[year = k]) + e_it
#   a_i = country FE, g_t = year FE, SEs clustered by country.
#   b_k = (PEPFAR - control gap in year k) - (the same gap in 2003).
#
# Sections:
#   1. Data and samples (in analysis/econ_setup.R)
#   2. Main event studies (unadjusted / adjusted), pre-trend tests, average effects
#   3. Linear time trends
#   4. Staggered timing (Sun and Abraham)
#   5. Robustness (winsorized, levels, balanced, Crown window, no China/India)
#   6. Figures and tables
#
# Run from the repository root:  Rscript analysis/parallel_trends/event_study_economic.R
# Writes everything to tables/crown_event_studies/robustness/

# Data, samples, outcomes, covariates, run_es(), tidy_es(), avg_coefs(), pre_test()
# and the plot style all come from the shared setup file
source("analysis/econ_setup.R")

out_dir <- "tables/crown_event_studies/robustness"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 2. Main event studies
# =====================================================

# The plain TWFE DiD on the same sample, to compare with the event-study average
run_did <- function(data, y) {
  m <- feols(as.formula(paste0(y, " ~ treat:post | iso3c + year")),
             data = data, cluster = ~iso3c)
  c(est = unname(coef(m)["treat:post"]), se = unname(se(m)["treat:post"]))
}

summarise_es <- function(m, data, y, group, spec) {
  pre  <- pre_test(m)
  apre <- avg_coefs(m, pre_pat)
  apst <- avg_coefs(m, post_pat)
  tibble(outcome = y, group = group, spec = spec,
         n_obs = nobs(m), n_countries = unname(m$fixef_sizes["iso3c"]),
         pre_wald_p = pre["p"], avg_pre = apre["est"],
         avg_post = apst["est"], avg_post_se = apst["se"])
}

es_coefs <- list()  # year-by-year coefficients, every spec
es_summ  <- list()  # one row per model

for (g in groups) for (y in names(outcomes)) for (spec in c("Unadjusted", "Adjusted")) {
  d <- samples[[g]]
  m <- run_es(d, y, adjusted = spec == "Adjusted")
  es_coefs[[length(es_coefs) + 1]] <- tidy_es(m) |>
    mutate(outcome = y, group = g, spec = spec)
  s <- summarise_es(m, d, y, g, spec)
  if (spec == "Unadjusted") {
    did <- run_did(d, y)
    s <- s |> mutate(did_twfe = did["est"], did_twfe_se = did["se"])
  }
  es_summ[[length(es_summ) + 1]] <- s
}

# =====================================================
# 3. Linear time trends (unadjusted)
# =====================================================
# (a) Slide version: Y = a_i + g_t + theta*(PEPFAR x t) + delta*(PEPFAR x Post).
#     theta = PEPFAR's extra drift per year over the whole period, delta = the jump
#     on top of it. If the effect builds gradually, theta can absorb it.
# (b) Country-specific linear trends: every country gets its own slope.
# (c) Pre-trend extrapolation: estimate theta on 1990-2003 only, subtract that line
#     from the PEPFAR countries' outcomes, then re-run the event study. Post
#     coefficients then read as "deviation from where the pre-2004 trend was heading".
#     (SEs here ignore the uncertainty in theta, so they are a bit too small.)

trend_rows <- list()
detr_coefs <- list()

for (g in groups) for (y in names(outcomes)) {
  d <- samples[[g]]

  ma <- feols(as.formula(paste0(y, " ~ treat:post + treat:t | iso3c + year")),
              data = d, cluster = ~iso3c)
  mb <- feols(as.formula(paste0(y, " ~ treat:post | iso3c[t] + year")),
              data = d, cluster = ~iso3c)
  mp <- feols(as.formula(paste0(y, " ~ treat:t | iso3c + year")),
              data = filter(d, year <= ref_year), cluster = ~iso3c)
  theta_pre <- unname(coef(mp)["treat:t"])

  d_dt <- d |> mutate(y_dt = .data[[y]] - theta_pre * treat * t)
  mc <- run_es(d_dt, "y_dt")
  detr_coefs[[length(detr_coefs) + 1]] <- tidy_es(mc) |>
    mutate(outcome = y, group = g, spec = "Pre-trend removed")
  ac <- avg_coefs(mc, post_pat)

  trend_rows[[length(trend_rows) + 1]] <- tibble(
    outcome = y, group = g,
    slide_delta = coef(ma)["treat:post"], slide_delta_se = se(ma)["treat:post"],
    slide_theta = coef(ma)["treat:t"],    slide_theta_se = se(ma)["treat:t"],
    ctrytrend_delta = coef(mb)["treat:post"], ctrytrend_delta_se = se(mb)["treat:post"],
    pre_theta = theta_pre, pre_theta_se = unname(se(mp)["treat:t"]),
    detrended_avg_post = ac["est"], detrended_avg_post_se = ac["se"]
  )
}

# =====================================================
# 4. Staggered timing: Sun and Abraham (2021)
# =====================================================
# Crown treats every PEPFAR country as starting in 2004. In the data, 59 of 90
# started in 2004 and the rest between 2005 and 2016. Here event time is counted from
# each country's own first disbursement; controls are never treated. Sun-Abraham
# avoids the bias plain TWFE has when timing is staggered.

# fixest reports Sun-Abraham coefficients averaged across cohorts for each event time,
# but its vcov() is for the underlying cohort x event-time coefficients. To test the
# leads jointly we rebuild the averaging (weights = each cohort's share of treated
# observations at that event time, as fixest does) and carry the vcov through it.
sa_lead_test <- function(m, d, leads = -14:-2) {
  b <- m$coefficients; nm <- names(b)
  per <- as.integer(sub("year::(-?\\d+):cohort::\\d+", "\\1", nm))
  coh <- as.integer(sub(".*cohort::", "", nm))
  cnt <- d[obs(m), ] |> filter(cohort < 10000) |> mutate(rel = year - cohort) |>
    count(cohort, rel)
  n <- cnt$n[match(paste(coh, per), paste(cnt$cohort, cnt$rel))]
  leads <- intersect(leads, per)
  W <- t(sapply(leads, function(p) { w <- ifelse(per == p, n, 0); w / sum(w) }))
  ba <- drop(W %*% b); Va <- W %*% vcov(m) %*% t(W)
  q <- length(leads); G <- unname(m$fixef_sizes["iso3c"])
  f_stat <- drop(t(ba) %*% solve(Va) %*% ba) / q
  pf(f_stat, q, G - 1, lower.tail = FALSE)
}

sa_coefs <- list()
sa_rows  <- list()

for (g in groups) for (y in names(outcomes)) {
  d <- samples[[g]] |>
    mutate(cohort = ifelse(treat == 1, pepfar_first_year, 10000))  # 10000 = never
  m <- feols(as.formula(paste0(y, " ~ sunab(cohort, year) | iso3c + year")),
             data = d, cluster = ~iso3c)
  ct <- coeftable(m)  # aggregated by event time
  keep <- grepl("^year::-?\\d+$", rownames(ct))
  sa_coefs[[length(sa_coefs) + 1]] <- tibble(
      event_time = as.integer(sub("year::", "", rownames(ct)[keep])),
      est = ct[keep, "Estimate"], se = ct[keep, "Std. Error"]) |>
    bind_rows(tibble(event_time = -1L, est = 0, se = 0)) |>
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se,
           outcome = y, group = g, spec = "Sun-Abraham (own start year)")
  att <- coeftable(summary(m, agg = "ATT"))["ATT", ]
  sa_rows[[length(sa_rows) + 1]] <- tibble(
    outcome = y, group = g, spec = "Sun-Abraham", n_obs = nobs(m),
    n_countries = unname(m$fixef_sizes["iso3c"]), pre_wald_p = sa_lead_test(m, d),
    avg_post = att[["Estimate"]], avg_post_se = att[["Std. Error"]])
}

# =====================================================
# 5. Robustness (unadjusted event studies, summarised)
# =====================================================

winsor <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

# Countries observed in at least half the pre years and half the post years for
# this outcome, so the comparison is not driven by countries entering/leaving
balanced <- function(d, y) {
  keep <- d |>
    filter(!is.na(.data[[y]])) |>
    group_by(iso3c) |>
    summarise(n_pre = sum(year < 2004), n_post = sum(year >= 2004), .groups = "drop") |>
    filter(n_pre >= 7, n_post >= 10) |>   # half of 1990-2003 (14) and 2004-2022 (19)
    pull(iso3c)
  filter(d, iso3c %in% keep)
}

rob_rows <- list()
add_rob <- function(d, y, g, spec, y_label = y) {
  m <- run_es(d, y)
  rob_rows[[length(rob_rows) + 1]] <<- summarise_es(m, d, y_label, g, spec)
}

growth_vars <- c("gdp_pc_growth", "gdp_growth")
level_map   <- c(log_gdp_pc = "gdp_pc", log_gdp = "gdp_bn")

for (g in groups) {
  crown <- make_sample(g, last_year = 2018, control = "grp_control")
  no_ci <- samples[[g]] |> filter(!iso3c %in% c("CHN", "IND"))
  wins  <- samples[[g]] |> mutate(across(all_of(growth_vars), winsor))
  for (y in names(outcomes)) {
    add_rob(crown, y, g, "Crown window 1990-2018")
    add_rob(no_ci, y, g, "Drop China & India")
    add_rob(balanced(samples[[g]], y), y, g, "Balanced countries")
    if (y %in% growth_vars) add_rob(wins, y, g, "Winsorized 1/99")
    if (y %in% names(level_map)) add_rob(samples[[g]], level_map[[y]], g,
                                         "Levels (2022 USD)", y_label = y)
  }
}

# =====================================================
# 6. Tables and figures
# =====================================================

es_coefs   <- bind_rows(es_coefs, detr_coefs)
es_summary <- bind_rows(es_summ, sa_rows, rob_rows) |>
  mutate(outcome_label = outcomes[outcome])
trends     <- bind_rows(trend_rows) |> mutate(outcome_label = outcomes[outcome])
sa_coefs   <- bind_rows(sa_coefs)

write_csv(es_coefs,   file.path(out_dir, "es_coefficients.csv"))
write_csv(sa_coefs,   file.path(out_dir, "sunab_coefficients.csv"))
write_csv(es_summary, file.path(out_dir, "es_summary.csv"))
write_csv(trends,     file.path(out_dir, "linear_trends.csv"))

# (1) One figure per outcome: All PEPFAR | COP, unadjusted vs adjusted
main_coefs <- es_coefs |> filter(spec %in% c("Unadjusted", "Adjusted")) |>
  mutate(spec = factor(spec, c("Unadjusted", "Adjusted")))

for (y in names(outcomes)) {
  ps <- lapply(groups, function(g) {
    pv <- es_summary |> filter(outcome == y, group == g,
                               spec %in% c("Unadjusted", "Adjusted"))
    sub <- paste0("Pre-2004 joint test p: unadjusted ",
                  fmt_p(pv$pre_wald_p[pv$spec == "Unadjusted"]),
                  ", adjusted ", fmt_p(pv$pre_wald_p[pv$spec == "Adjusted"]))
    ggplot(filter(main_coefs, outcome == y, group == g),
           aes(year, est, color = spec)) +
      es_layer(2003.5) +
      scale_color_manual(values = pal[1:2]) +
      labs(title = paste(g, "vs. controls"), subtitle = sub,
           x = NULL, y = "Gap vs. 2003") +
      theme_es
  })
  p <- (ps[[1]] | ps[[2]]) + plot_layout(guides = "collect") +
    plot_annotation(
      title = outcomes[[y]],
      caption = paste("Coefficients on PEPFAR x year, 2003 = 0, with 95% CIs.",
                      "Country and year FE; SEs clustered by country.",
                      "Adjusted adds baseline covariates x year."),
      theme = theme(plot.title = element_text(face = "bold", size = 13, color = ink),
                    plot.caption = element_text(color = ink2, hjust = 0),
                    legend.position = "top")) &
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, paste0("es_", y, ".png")), p,
         width = 11, height = 4.2, dpi = 300, bg = "white")
}

# (2) Overview grid: all 8 outcomes, unadjusted, All PEPFAR vs COP
lab_fac <- function(x) factor(outcomes[x], levels = outcomes)
p_grid <- main_coefs |> filter(spec == "Unadjusted") |>
  mutate(outcome = lab_fac(outcome)) |>
  ggplot(aes(year, est, color = group)) +
  es_layer(2003.5) +
  facet_wrap(~outcome, scales = "free_y", ncol = 2) +
  scale_color_manual(values = pal[1:2]) +
  labs(title = "Event studies, Crown economic outcomes (unadjusted)",
       subtitle = "PEPFAR x year coefficients relative to 2003, 95% CIs; controls = grp_control_kates",
       x = NULL, y = "Gap vs. 2003") +
  theme_es
ggsave(file.path(fig_dir, "es_overview_unadjusted.png"), p_grid,
       width = 10, height = 11, dpi = 300, bg = "white")

# (3) Common 2004 start vs. Sun-Abraham own start year, by event time
cmp_sa <- bind_rows(
  main_coefs |> filter(spec == "Unadjusted") |>
    transmute(event_time = year - 2004, est, lo, hi, outcome, group,
              spec = "Common 2004 start (TWFE)"),
  sa_coefs |> select(event_time, est, lo, hi, outcome, group, spec)) |>
  filter(event_time >= -14, event_time <= 18) |>
  mutate(outcome = lab_fac(outcome))

for (g in groups) {
  p <- cmp_sa |> filter(group == g) |>
    ggplot(aes(event_time, est, color = spec)) +
    es_layer(-0.5) +
    facet_wrap(~outcome, scales = "free_y", ncol = 2) +
    scale_color_manual(values = pal[c(1, 3)]) +
    labs(title = paste0("Staggered timing check: ", g, " vs. controls"),
         subtitle = "Event time 0 = first PEPFAR year (2004 for the common-start version); -1 = 0",
         x = "Years since PEPFAR start", y = "Gap vs. year before start") +
    theme_es
  ggsave(file.path(fig_dir, paste0("es_sunab_", gsub(" ", "_", tolower(g)), ".png")),
         p, width = 10, height = 11, dpi = 300, bg = "white")
}

# (4) Raw vs. pre-trend removed
for (g in groups) {
  p <- es_coefs |> filter(group == g, spec %in% c("Unadjusted", "Pre-trend removed")) |>
    mutate(outcome = lab_fac(outcome),
           spec = factor(spec, c("Unadjusted", "Pre-trend removed"))) |>
    ggplot(aes(year, est, color = spec)) +
    es_layer(2003.5) +
    facet_wrap(~outcome, scales = "free_y", ncol = 2) +
    scale_color_manual(values = pal[c(1, 2)]) +
    labs(title = paste0("Linear pre-trend check: ", g, " vs. controls"),
         subtitle = "Pre-trend removed: PEPFAR's 1990-2003 linear trend vs. controls extrapolated forward and subtracted",
         x = NULL, y = "Gap vs. 2003") +
    theme_es
  ggsave(file.path(fig_dir, paste0("es_detrended_", gsub(" ", "_", tolower(g)), ".png")),
         p, width = 10, height = 11, dpi = 300, bg = "white")
}

# (5) Average post-2004 effect across every specification
spec_order <- c("Unadjusted", "Adjusted", "Sun-Abraham", "Crown window 1990-2018",
                "Drop China & India", "Balanced countries", "Winsorized 1/99",
                "Levels (2022 USD)")
forest <- es_summary |>
  filter(spec != "Levels (2022 USD)") |>   # different units from the log version
  mutate(spec = factor(spec, rev(spec_order)), outcome = lab_fac(outcome),
         lo = avg_post - 1.96 * avg_post_se, hi = avg_post + 1.96 * avg_post_se)
p_forest <- ggplot(forest, aes(avg_post, spec, color = group)) +
  geom_vline(xintercept = 0, color = ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.5,
                orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  facet_wrap(~outcome, scales = "free", ncol = 2) +  # free y drops rows a check skips
  scale_color_manual(values = pal[1:2]) +
  labs(title = "Average post-2004 effect across specifications",
       subtitle = "Mean of the 2004-2022 event-study coefficients (Sun-Abraham: its ATT), 95% CIs",
       x = "Average post-period gap vs. 2003", y = NULL) +
  theme_es + theme(panel.grid.major.x = element_line(color = grid_col, linewidth = 0.3),
                   panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "robustness_avg_post.png"), p_forest,
       width = 10, height = 11, dpi = 300, bg = "white")

cat("Done. Results in", out_dir, "\n")

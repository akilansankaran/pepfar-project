# Adjusted event studies for the Crown et al. (2023) economic outcomes:
# coefficient tables (LaTeX) and graphs (PDF for LaTeX, PNG for quick viewing)
# Author: Anthony Babu
#
# Model (one per outcome x treatment group):
#   Y_it = a_i + g_t + sum_{k != 2003} b_k (PEPFAR_i x 1[year = k]) + X_i' d_t + e_it
#   a_i = country FE, g_t = year FE, X_i' d_t = each of the 11 baseline covariates with
#   its own coefficient in every year (income group as income group x year FE).
#   b_k = (PEPFAR - control gap in year k) - (the same gap in 2003). SEs clustered by country.
# This is the year-by-year version of Panel B (TWFE) in did_econ.R.
#
# Run from the repository root:  Rscript analysis/parallel_trends/event_study_econ_adjusted.R
# Writes to tables/crown_event_studies/adjusted/{tex, figures, csv}/
# The LaTeX files refer to each other and to the figures by their path from the
# tables/ folder (e.g. crown_event_studies/adjusted/figures/...), so they work both
# from tables/crown_results.tex and in Overleaf if the same folders are uploaded.
# Labels start with crown_ so they never clash with the mortality tables.

source("analysis/econ_setup.R")

tex_root <- "crown_event_studies/adjusted"   # path from tables/ (and from the Overleaf root)
out_dir <- file.path("tables", tex_root)
tex_dir <- file.path(out_dir, "tex")
fig_dir <- file.path(out_dir, "figures")
csv_dir <- file.path(out_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Estimate
# =====================================================

coefs <- list()  # year-by-year coefficients
summ  <- list()  # one row per model: N, countries, pre-trend test, average post effect

for (g in groups) for (y in names(outcomes)) {
  d <- samples[[g]]
  m <- run_es(d, y, adjusted = TRUE)
  coefs[[length(coefs) + 1]] <- tidy_es(m) |> mutate(outcome = y, group = g)
  pre <- pre_test(m)
  ap <- avg_coefs(m, post_pat)
  G  <- n_countries(m, d)
  summ[[length(summ) + 1]] <- tibble(
    outcome = y, group = g, n_obs = nobs(m), n_countries = G,
    pre_p = pre[["p"]], avg_post = ap[["est"]], avg_post_se = ap[["se"]],
    avg_post_p = 2 * pt(-abs(ap[["est"]] / ap[["se"]]), df = G - 1))
}
coefs <- bind_rows(coefs)
summ  <- bind_rows(summ)
write_csv(coefs, file.path(csv_dir, "es_adjusted_coefficients.csv"))
write_csv(summ,  file.path(csv_dir, "es_adjusted_summary.csv"))

# =====================================================
# 2. LaTeX tables
# =====================================================

ys        <- names(outcomes)
pre_yrs   <- 1990:2002
post_yrs  <- 2004:2022

# Coefficient or SE for one outcome/group/year; blank if that year has no data
es_cell <- function(g, y, yr, what) {
  r <- filter(coefs, group == g, outcome == y, year == yr)
  if (nrow(r) == 0) return("")
  if (what == "coef") tex_coef(r$est, r$p) else tex_se(r$se)
}
s_cell <- function(g, y, what) {
  r <- filter(summ, group == g, outcome == y)
  switch(what,
         pre_p = if (r$pre_p < 0.001) "$<$0.001" else sprintf("%.3f", r$pre_p),
         avg = tex_coef(r$avg_post, r$avg_post_p), avg_se = tex_se(r$avg_post_se),
         n = tex_n(r$n_obs), countries = tex_n(r$n_countries))
}

es_notes <- paste(
  "Coefficients on treatment group $\\times$ year, relative to 2003 (omitted);",
  "standard errors clustered by country in parentheses. All models include country",
  "and year fixed effects, the 11 baseline (2003) covariates of Crown et al.\\ (2023)",
  "each with its own coefficient in every year, and World Bank income group",
  "$\\times$ year fixed effects. Sample 1990--2022; controls are the 64 comparison",
  "countries valid through 2022 (\\texttt{grp\\_control\\_kates}). Blank cells:",
  "no data for that outcome in that year. Disengagement is in percentage points.",
  "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.")

# ---- Separate tables: one group, columns = 8 outcomes, split pre / post ----
# (66 coefficient rows don't fit on one page, so the pre-2004 and post-2004 years
#  are two tables; the pre table ends with the parallel-trends test, the post table
#  with the average effect.)
for (g in groups) {
  year_rows <- function(yrs) {
    unlist(lapply(yrs, function(yr) list(
      c(as.character(yr), sapply(ys, \(y) es_cell(g, y, yr, "coef"))),
      c("",               sapply(ys, \(y) es_cell(g, y, yr, "se"))))), recursive = FALSE)
  }
  foot <- list(c("Observations", sapply(ys, \(y) s_cell(g, y, "n"))),
               c("Countries",    sapply(ys, \(y) s_cell(g, y, "countries"))))
  head <- list(c("Year", sprintf("(%d)", seq_along(ys))), c("", outcomes_tex[ys]))
  align <- paste0("l", strrep("c", length(ys)))
  gf <- group_file[[g]]

  write_tex_table(
    file    = file.path(tex_dir, paste0("es_adjusted_", gf, "_pre.tex")),
    caption = paste0("Adjusted event study, ", g, " vs.\\ controls: pre-period (1990--2002)"),
    label   = paste0("tab:crown_es_", gf, "_pre"),
    head    = head,
    rows    = c(year_rows(pre_yrs), list("\\midrule",
                c("Pre-2004 joint test ($p$)", sapply(ys, \(y) s_cell(g, y, "pre_p")))), foot),
    notes   = paste(es_notes, "The joint test asks whether all 1990--2002",
                    "coefficients are zero; a small $p$ means trends were not parallel",
                    "before PEPFAR."),
    align   = align, size = "\\small")

  write_tex_table(
    file    = file.path(tex_dir, paste0("es_adjusted_", gf, "_post.tex")),
    caption = paste0("Adjusted event study, ", g, " vs.\\ controls: post-period (2004--2022)"),
    label   = paste0("tab:crown_es_", gf, "_post"),
    head    = head,
    rows    = c(year_rows(post_yrs), list("\\midrule",
                c("Average 2004--2022", sapply(ys, \(y) s_cell(g, y, "avg"))),
                c("",                   sapply(ys, \(y) s_cell(g, y, "avg_se")))), foot),
    notes   = paste(es_notes, "Average 2004--2022 is the mean of the post-period",
                    "coefficients, with its standard error."),
    align   = align, size = "\\small")
}

# ---- Combined tables: one per outcome, All PEPFAR and COP side by side ----
combined_files <- character(0)
for (y in ys) {
  rows <- lapply(c(pre_yrs, post_yrs), function(yr) {
    c(as.character(yr), es_cell("All PEPFAR", y, yr, "coef"), es_cell("All PEPFAR", y, yr, "se"),
      es_cell("COP", y, yr, "coef"), es_cell("COP", y, yr, "se"))
  })
  # put the reference year in its place, between 2002 and 2004
  rows <- c(rows[seq_along(pre_yrs)], list(c("2003", "\\multicolumn{2}{c}{(reference)}",
                                             "\\multicolumn{2}{c}{(reference)}")),
            rows[-seq_along(pre_yrs)])
  # drop years with no data for either group (e.g. employment starts in 1991)
  rows <- Filter(\(r) any(r[-1] != ""), rows)
  sc <- function(what) c(s_cell("All PEPFAR", y, what), s_cell("COP", y, what))
  two <- function(label, a, b) c(label, sprintf("\\multicolumn{2}{c}{%s}", a),
                                 sprintf("\\multicolumn{2}{c}{%s}", b))
  rows <- c(rows, list("\\midrule",
    two("Pre-2004 joint test ($p$)", sc("pre_p")[1], sc("pre_p")[2]),
    c("Average 2004--2022", s_cell("All PEPFAR", y, "avg"), s_cell("All PEPFAR", y, "avg_se"),
      s_cell("COP", y, "avg"), s_cell("COP", y, "avg_se")),
    two("Observations", sc("n")[1], sc("n")[2]),
    two("Countries", sc("countries")[1], sc("countries")[2])))
  f <- paste0("es_adjusted_combined_", y, ".tex")
  combined_files <- c(combined_files, f)
  write_tex_table(
    file    = file.path(tex_dir, f),
    caption = paste0("Adjusted event study: ", outcomes_tex[[y]],
                     ", All PEPFAR and COP"),
    label   = paste0("tab:crown_es_combined_", y),
    head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR vs.\\ controls}",
                     "\\multicolumn{2}{c}{COP vs.\\ controls}"),
                   "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
                   c("Year", "Coef.", "SE", "Coef.", "SE")),
    rows    = rows,
    notes   = es_notes,
    align   = "lcccc", resize = FALSE, size = "\\small")
}
# One file that pulls in all 8 combined tables, each on its own page
writeLines(paste0("\\input{", tex_root, "/tex/", sub("\\.tex$", "", combined_files), "}\n\\clearpage"),
           file.path(tex_dir, "es_adjusted_combined_all.tex"))

# =====================================================
# 3. Graphs
# =====================================================
# Color follows the group everywhere: All PEPFAR = blue, COP = orange.

group_col <- c("All PEPFAR" = pal[1], "COP" = pal[2])
cap <- paste0("Coefficients on treatment x year, 2003 = 0, with 95% CIs. Country and year FE,\n",
              "baseline covariates x year, income group x year FE; SEs clustered by country.")
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, bg = "white")
}
sub_for <- function(y, gs) {
  paste(sapply(gs, function(g) {
    r <- filter(summ, outcome == y, group == g)
    p_txt <- if (r$pre_p < 0.001) "p < 0.001" else sprintf("p = %.3f", r$pre_p)
    sprintf("%s: pre-2004 test %s, avg. post = %.2f (%.2f)",
            g, p_txt, r$avg_post, r$avg_post_se)
  }), collapse = "\n")
}

plot_es <- function(dat, title, subtitle, legend = TRUE) {
  ggplot(dat, aes(year, est, color = group)) +
    es_layer(2003.5) +
    scale_color_manual(values = group_col) +
    labs(title = title, subtitle = subtitle, caption = cap,
         x = NULL, y = "Gap vs. 2003") +
    theme_es +
    theme(plot.caption = element_text(color = ink2, hjust = 0, size = 7),
          legend.position = if (legend) "top" else "none")
}

# (a) One graph per outcome, both groups together
# (b) One graph per outcome per group (separate)
for (y in ys) {
  save_both(plot_es(filter(coefs, outcome == y), outcomes[[y]], sub_for(y, groups)),
            paste0("es_adj_", y, "_combined"), 7, 4.5)
  for (g in groups) {
    save_both(plot_es(filter(coefs, outcome == y, group == g),
                      paste0(outcomes[[y]], ": ", g, " vs. controls"),
                      sub_for(y, g), legend = FALSE),
              paste0("es_adj_", y, "_", group_file[[g]]), 7, 4.2)
  }
}

# (c) All 8 outcomes on one page: both groups, and each group alone
grid_plot <- function(dat, title, legend = TRUE) {
  dat |> mutate(outcome = factor(outcomes[outcome], levels = outcomes)) |>
    ggplot(aes(year, est, color = group)) +
    es_layer(2003.5) +
    facet_wrap(~outcome, scales = "free_y", ncol = 2) +
    scale_color_manual(values = group_col) +
    labs(title = title, caption = cap, x = NULL, y = "Gap vs. 2003") +
    theme_es +
    theme(plot.caption = element_text(color = ink2, hjust = 0, size = 7),
          legend.position = if (legend) "top" else "none")
}
save_both(grid_plot(coefs, "Adjusted event studies: All PEPFAR and COP vs. controls"),
          "es_adj_grid_combined", 7.5, 9.5)
for (g in groups) {
  save_both(grid_plot(filter(coefs, group == g),
                      paste0("Adjusted event studies: ", g, " vs. controls"), legend = FALSE),
            paste0("es_adj_grid_", group_file[[g]]), 7.5, 9.5)
}

# =====================================================
# 4. LaTeX figure blocks
# =====================================================
fig_env <- function(file, caption, label, width = "\\textwidth") {
  c("\\begin{figure}[htbp]", "\\centering",
    sprintf("\\includegraphics[width=%s]{%s/figures/%s.pdf}", width, tex_root, file),
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    "\\end{figure}")
}
tex_lab <- function(y) gsub("%", "\\\\%", outcomes[[y]])

# Combined: grid first, then one per outcome (two per page)
lines <- c(fig_env("es_adj_grid_combined",
                   "Adjusted event studies, all outcomes: All PEPFAR (blue) and COP (orange)",
                   "fig:crown_es_grid_combined"), "\\clearpage")
for (i in seq_along(ys)) {
  y <- ys[i]
  lines <- c(lines, fig_env(paste0("es_adj_", y, "_combined"),
                            paste0("Adjusted event study: ", tex_lab(y)),
                            paste0("fig:crown_es_combined_", y)))
  if (i %% 2 == 0) lines <- c(lines, "\\clearpage")
}
writeLines(lines, file.path(tex_dir, "figures_combined.tex"))

# Separate: per group, grid then one per outcome
for (g in groups) {
  gf <- group_file[[g]]
  lines <- c(fig_env(paste0("es_adj_grid_", gf),
                     paste0("Adjusted event studies, all outcomes: ", g, " vs.\\ controls"),
                     paste0("fig:crown_es_grid_", gf)), "\\clearpage")
  for (i in seq_along(ys)) {
    y <- ys[i]
    lines <- c(lines, fig_env(paste0("es_adj_", y, "_", gf),
                              paste0("Adjusted event study: ", tex_lab(y), ", ", g),
                              paste0("fig:crown_es_", gf, "_", y)))
    if (i %% 2 == 0) lines <- c(lines, "\\clearpage")
  }
  writeLines(lines, file.path(tex_dir, paste0("figures_", gf, ".tex")))
}

cat("Adjusted event studies done. Output in", out_dir, "\n")

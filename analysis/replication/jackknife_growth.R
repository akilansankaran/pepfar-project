# Jackknife: does any single country drive the panel growth results?
# Author: Anthony Babu
#
# The Equatorial Guinea finding (one control country carried the triple interaction)
# raises the general question: is any other result driven by one country? This script
# re-runs every regression in the Crown panel growth tables (Tables 3-4,
# twfe_growth_robustness.R) once for each country, leaving that country out, and records
# the coefficient and p-value each time ("leave-one-out" or jackknife).
#
#   Regressions: GDP per capita growth and GDP growth;
#     PEPFAR dummy x Post (run_twfe_post) and HIV x PEPFAR x Post (run_twfe_triple_post);
#     All PEPFAR and COP samples; unadjusted and adjusted (11 baseline covariates x Post);
#     raw, winsorized, winsorized log and winsorized IHS outcomes (as in Tables 3-4).
#   For each regression, every country in its estimation sample is dropped in turn, so
#   the number of runs equals the number of countries in that regression.
#
# Reported:
#   jackknife_growth.csv              every leave-one-out run (coefficient, SE, p)
#   jackknife_growth_summary.csv/.tex share of runs significant at 10%, and the range of
#                                     the coefficient with the country whose removal moves
#                                     it most in each direction
#   jackknife_matrix_*.pdf/.png       the full matrix: one row per dropped country, one
#                                     column per regression (raw outcome), coloured by
#                                     significance
#
# Run from the repository root:  Rscript analysis/replication/jackknife_growth.R
# Writes to tables/crown_panel/jackknife/

source("analysis/econ_setup.R")

out_dir <- "tables/crown_panel/jackknife"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

growth_outcomes <- c(gdp_pc_growth = "GDP per capita growth", gdp_growth = "GDP growth")
transforms <- c(raw = "Raw", win = "Winsorized", winlog = "Winsorized, log", winihs = "Winsorized, IHS")
specs      <- c(unadj = "Unadj.", adj = "Adj.")
treat_lab  <- c(treat = "PEPFAR $\\times$ Post", triple = "HIV $\\times$ PEPFAR $\\times$ Post")

# Same transformations as twfe_growth_robustness.R (cut-offs from the full panel)
panel_t <- panel
for (y in names(growth_outcomes)) {
  w <- winsor(panel[[y]])
  panel_t[[paste0(y, "__raw")]]    <- panel[[y]]
  panel_t[[paste0(y, "__win")]]    <- w
  panel_t[[paste0(y, "__winlog")]] <- 100 * log1p(w / 100)
  panel_t[[paste0(y, "__winihs")]] <- asinh(w)
}

fit <- function(d, yvar, tv, adjusted) {
  if (tv == "treat") {
    m <- run_twfe_post(d, yvar, tvar = "treat", adjusted = adjusted)
    ct <- coeftable(m)["treat:post", ]
  } else {
    m <- run_twfe_triple_post(d, yvar, adjusted = adjusted)
    ct <- coeftable(m)["hiv_treat_post", ]
  }
  list(m = m, est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]])
}

# ---- Full-sample estimates and every leave-one-out run ----
full <- list(); runs <- list()
for (g in groups) {
  d0 <- make_sample(g, data = panel_t) |> filter(!is.na(bl_hiv_prev))
  for (y in names(growth_outcomes)) for (tr in names(transforms))
    for (tv in names(treat_lab)) for (sp in names(specs)) {
      yvar <- paste0(y, "__", tr)
      f0 <- fit(d0, yvar, tv, sp == "adj")
      key <- tibble(outcome = y, transform = tr, treatment = tv, group = g, spec = sp)
      full[[length(full) + 1]] <- bind_cols(key, tibble(est = f0$est, se = f0$se, p = f0$p))
      in_sample <- sort(unique(d0$iso3c[obs(f0$m)]))   # only countries this regression uses
      for (ctry in in_sample) {
        f <- fit(filter(d0, iso3c != ctry), yvar, tv, sp == "adj")
        runs[[length(runs) + 1]] <- bind_cols(key, tibble(dropped = ctry, est = f$est, se = f$se, p = f$p))
      }
    }
}
full <- bind_rows(full); runs <- bind_rows(runs)
write_csv(runs, file.path(out_dir, "jackknife_growth.csv"))

summ <- runs |>
  group_by(outcome, transform, treatment, group, spec) |>
  summarise(n_runs = n(), share_sig10 = mean(p < 0.10),
            min_est = min(est), min_drop = dropped[which.min(est)],
            max_est = max(est), max_drop = dropped[which.max(est)],
            not_sig_when_dropped = paste(dropped[p >= 0.10], collapse = " "), .groups = "drop") |>
  left_join(full, by = c("outcome", "transform", "treatment", "group", "spec"))
write_csv(summ, file.path(out_dir, "jackknife_growth_summary.csv"))

# ---- LaTeX summary: rows = regression (treatment x group x spec), columns = transformation ----
pct <- function(x) sprintf("%.0f\\%%", 100 * x)
rows <- list()
for (k in seq_along(growth_outcomes)) {
  y <- names(growth_outcomes)[k]
  if (k > 1) rows[[length(rows) + 1]] <- "\\midrule"
  rows[[length(rows) + 1]] <- sprintf("\\multicolumn{7}{l}{\\textbf{Panel %s: %s}} \\\\",
                                      LETTERS[k], growth_outcomes[[y]])
  for (tv in names(treat_lab)) for (g in groups) for (sp in names(specs)) {
    r <- filter(summ, outcome == y, treatment == tv, group == g, spec == sp)
    raw <- filter(r, transform == "raw")
    rows[[length(rows) + 1]] <- c(
      sprintf("%s, %s, %s", treat_lab[[tv]], g, specs[[sp]]),
      tex_coef(raw$est, raw$p),
      sapply(names(transforms), \(tr) pct(r$share_sig10[r$transform == tr])),
      sprintf("%s to %s", formatC(raw$min_est, format = "f", digits = 2),
              formatC(raw$max_est, format = "f", digits = 2)))
  }
}
write_tex_table(
  file    = file.path(out_dir, "jackknife_growth_summary.tex"),
  caption = "Jackknife: share of leave-one-country-out regressions significant at 10\\%",
  label   = "tab:crown_jackknife",
  head    = list(c("", "Full sample", "\\multicolumn{4}{c}{Share significant at 10\\%}", "Range"),
                 "\\cmidrule(lr){3-6}",
                 c("Regression", "(raw)", "Raw", "Wins.", "Wins., log", "Wins., IHS", "(raw)")),
  rows    = rows,
  notes   = paste(
    "Each regression of Tables 3--4 is re-estimated once for every country in its sample, leaving",
    "that country out (about 150 runs for All PEPFAR, 90 for COP, fewer when adjusted). Full",
    "sample: the coefficient in Tables 3--4 with the raw outcome. Share significant at 10\\%:",
    "fraction of the leave-one-out runs in which the coefficient has $p<0.10$, for each outcome",
    "transformation. Range: smallest and largest coefficient across the leave-one-out runs (raw",
    "outcome). PEPFAR $\\times$ Post: two-way fixed effects DiD with the 0/1 PEPFAR (or COP) dummy.",
    "HIV $\\times$ PEPFAR $\\times$ Post: the triple interaction with 2003 HIV prevalence, which",
    "also includes HIV $\\times$ Post and PEPFAR $\\times$ Post. Adj.: the 11 baseline (2003)",
    "covariates of Crown et al.\\ (2023) interacted with Post. Country and year fixed effects;",
    "standard errors clustered by country. Countries without 2003 HIV prevalence (3 controls) are",
    "left out of every regression so that both treatments use the same countries."),
  align = "lcccccc", resize = TRUE)

# ---- The full matrix as a figure: one row per dropped country, one column per regression ----
col_order <- unlist(lapply(names(treat_lab), \(tv) lapply(groups, \(g) lapply(names(specs), \(sp)
  paste0(c(treat = "PEPFAR", triple = "Triple")[[tv]], "\n", g, "\n", specs[[sp]])))))
sig_class <- function(p) cut(p, c(-Inf, 0.01, 0.05, 0.10, Inf),
                             labels = c("p < 0.01", "p < 0.05", "p < 0.10", "not significant"))
for (y in names(growth_outcomes)) {
  m <- runs |>
    filter(outcome == y, transform == "raw") |>
    mutate(column = factor(paste0(c(treat = "PEPFAR", triple = "Triple")[treatment], "\n", group,
                                  "\n", specs[spec]), levels = col_order),
           sig = sig_class(p),
           label = formatC(est, format = "f", digits = 1))
  ctry_order <- sort(unique(m$dropped), decreasing = TRUE)
  m$dropped <- factor(m$dropped, levels = ctry_order)
  p <- ggplot(m, aes(column, dropped, fill = sig)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = 1.9) +
    scale_fill_manual(values = c("p < 0.01" = "#1f6fb2", "p < 0.05" = "#5fa2d9",
                                 "p < 0.10" = "#b7d6ef", "not significant" = "#f2c4b5"),
                      drop = FALSE, name = NULL) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = "Country left out",
         title = paste0("Jackknife matrix: ", growth_outcomes[[y]], " (raw)"),
         subtitle = "Each cell: coefficient when that country is dropped; colour: significance. Blank: country not in that regression's sample.") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(size = 5.5), plot.title.position = "plot")
  for (ext in c("pdf", "png"))
    ggsave(file.path(out_dir, paste0("jackknife_matrix_", y, ".", ext)), p,
           width = 7.5, height = 15, dpi = 200)
}

# LaTeX block for the two matrix figures
fig_tex <- unlist(lapply(names(growth_outcomes), \(y) c(
  "\\begin{figure}[p]", "\\centering",
  sprintf("\\includegraphics[height=0.92\\textheight]{crown_panel/jackknife/jackknife_matrix_%s.pdf}", y),
  sprintf("\\caption{Jackknife matrix, %s (raw outcome): coefficient with each country left out}",
          tolower(growth_outcomes[[y]])),
  sprintf("\\label{fig:crown_jackknife_%s}", y), "\\end{figure}", "")))
writeLines(fig_tex, file.path(out_dir, "jackknife_figures.tex"))

cat("Jackknife done:", nrow(runs), "leave-one-out regressions. Output in", out_dir, "\n")

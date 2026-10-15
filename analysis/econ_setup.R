# Shared setup for the economic-outcome DiD and event-study scripts
# Author: Anthony Babu
#
# Sourced by:
#   analysis/replication/did_econ.R
#   analysis/parallel_trends/event_study_econ_adjusted.R
#   analysis/parallel_trends/event_study_economic.R
# so the data, samples, outcomes and covariates are defined in exactly one place.

suppressPackageStartupMessages({
  library(dplyr)     # data editing
  library(readr)     # reading/writing CSVs
  library(tidyr)     # reshaping
  library(fixest)    # regressions with fixed effects
  library(ggplot2)   # plots
  library(patchwork) # putting plots side by side
})
setFixest_notes(FALSE)  # silence the "observations removed because of NA" messages

ref_year <- 2003  # the omitted year in the event studies; last pre-PEPFAR year

# =====================================================
# Data
# =====================================================

panel <- read_csv("data/intermediate/pepfar_panel_1990_2022_FIXED.csv",
                  show_col_types = FALSE) |>
  mutate(
    # Crown reports disengagement in percentage points, the panel stores a 0-1 share
    prim_disengage_girls = 100 * prim_disengage_girls,
    prim_disengage_boys  = 100 * prim_disengage_boys,
    # Levels of GDP are very skewed (China vs. Tuvalu), so the main version uses logs
    log_gdp_pc = log(gdp_pc),
    log_gdp    = log(gdp),
    gdp_bn     = gdp / 1e9,   # level in billions, used only as a robustness check
    post = as.integer(year >= 2004),  # Crown's "Time": 1 from 2004 on, every country
    t    = year - ref_year,           # linear time, 0 in 2003
    bl_wb_income = na_if(bl_wb_income, "")
  )

# Order: Crown's five original outcomes, then the team's three additions
outcomes <- c(
  gdp_pc_growth        = "GDP per capita growth (%)",
  emp_female           = "Female employment (% of 15+)",
  emp_male             = "Male employment (% of 15+)",
  prim_disengage_girls = "Girls' primary disengagement (pp)",
  prim_disengage_boys  = "Boys' primary disengagement (pp)",
  gdp_growth           = "GDP growth (%)",
  log_gdp_pc           = "log GDP per capita",
  log_gdp              = "log GDP"
)

# Short names for LaTeX table headers (no % signs, which LaTeX treats as comments)
outcomes_tex <- c(
  gdp_pc_growth        = "GDP pc growth",
  emp_female           = "Female emp.",
  emp_male             = "Male emp.",
  prim_disengage_girls = "Girls' diseng.",
  prim_disengage_boys  = "Boys' diseng.",
  gdp_growth           = "GDP growth",
  log_gdp_pc           = "log GDP pc",
  log_gdp              = "log GDP"
)

# The 11 baseline (2003) covariates from Crown Table 1, as built in the panel.
# Ten are numbers (or 0/1); income group is a category (L / LM / UM).
bl_cont <- c("bl_gdp_pc", "bl_pop_total", "bl_life_exp", "bl_fertility",
             "bl_urban_pct", "bl_sec_enroll_gross", "bl_hiv_prev",
             "bl_donor_health_nonpepfar_pc", "bl_dom_hexp_pc_ppp",
             "bl_us_hiv_pre2004")
bl_all <- c(bl_cont, "bl_wb_income")

# =====================================================
# Samples
# =====================================================
# Treatment group vs. controls. grp_control_kates is the control flag that is valid
# through 2022 (drops Venezuela, Colombia, Panama); grp_control is Crown's 2004-2018 one.
# Every sample gets one column, treat, that is 1 for the treatment group and 0 for controls.
make_sample <- function(group, last_year = 2022, control = "grp_control_kates",
                        data = panel) {
  tv <- c("All PEPFAR" = "grp_pepfar_any", "COP" = "grp_pepfar_cop")[[group]]
  data |>
    filter(year <= last_year, .data[[tv]] == 1 | .data[[control]] == 1) |>
    mutate(treat = .data[[tv]])
}

groups  <- c("All PEPFAR", "COP")
samples <- setNames(lapply(groups, make_sample), groups)
group_file <- c("All PEPFAR" = "all_pepfar", "COP" = "cop")  # for file names

# =====================================================
# Event-study helpers
# =====================================================

# Y_it = a_i + g_t + sum_{k != 2003} b_k (treat_i x 1[year = k]) + e_it
# Adjusted: each baseline covariate also gets its own coefficient in every year
# (year[x1, x2, ...] in fixest), and income group gets income-group x year FE.
run_es <- function(data, y, adjusted = FALSE) {
  fe <- if (adjusted) {
    paste0("iso3c + bl_wb_income^year + year[", paste(bl_cont, collapse = ", "), "]")
  } else {
    "iso3c + year"
  }
  f <- as.formula(paste0(y, " ~ i(year, treat, ref = ", ref_year, ") | ", fe))
  feols(f, data = data, cluster = ~iso3c)
}

# Pull the year-by-year coefficients into a tidy table (2003 added back as 0)
tidy_es <- function(m) {
  ct <- coeftable(m)
  keep <- grepl("^year::\\d+:treat$", rownames(ct))
  tibble(year = as.integer(sub("year::(\\d+):treat", "\\1", rownames(ct)[keep])),
         est  = ct[keep, "Estimate"],
         se   = ct[keep, "Std. Error"],
         p    = ct[keep, "Pr(>|t|)"]) |>
    bind_rows(tibble(year = ref_year, est = 0, se = 0, p = NA_real_)) |>
    arrange(year) |>
    mutate(lo = est - 1.96 * se, hi = est + 1.96 * se)
}

pre_pat  <- "^year::(199\\d|200[0-2]):treat$"        # 1990-2002
post_pat <- "^year::(200[4-9]|20[1-2]\\d):treat$"    # 2004-2029 (data stop in 2022)

# Average of a set of coefficients, with its standard error (w' V w)
avg_coefs <- function(m, pattern) {
  b <- coef(m); V <- vcov(m)
  idx <- grep(pattern, names(b))
  w <- rep(1 / length(idx), length(idx))
  c(est = sum(w * b[idx]), se = sqrt(drop(t(w) %*% V[idx, idx] %*% w)))
}

# Joint test that every pre-2004 coefficient is zero. A small p-value means the
# groups were already moving apart before PEPFAR, i.e. parallel trends looks violated.
pre_test <- function(m) {
  w <- wald(m, keep = pre_pat, print = FALSE)
  c(stat = unname(w$stat), p = unname(w$p))
}

n_countries <- function(m, data) length(unique(data$iso3c[obs(m)]))

# =====================================================
# Plot style
# =====================================================
# Colors: first three slots of a colorblind-validated categorical palette
pal <- c("#2a78d6", "#eb6834", "#1baf7a")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"

theme_es <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = grid_col, linewidth = 0.3),
        axis.text = element_text(color = ink2),
        axis.title = element_text(color = ink2),
        plot.title = element_text(color = ink, face = "bold"),
        plot.subtitle = element_text(color = ink2),
        strip.text = element_text(color = ink, face = "bold", hjust = 0),
        legend.position = "top", legend.justification = "left",
        legend.title = element_blank(), legend.text = element_text(color = ink))

# Zero line, dashed line between 2003 and 2004, CI whiskers, dots
es_layer <- function(x_break) {
  list(geom_hline(yintercept = 0, color = ink2, linewidth = 0.4),
       geom_vline(xintercept = x_break, color = ink2, linetype = "dashed",
                  linewidth = 0.4),
       geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5,
                     position = position_dodge(width = 0.6)),
       geom_point(size = 1.6, position = position_dodge(width = 0.6)))
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# =====================================================
# LaTeX table helpers
# =====================================================
# Stars: * p<0.10, ** p<0.05, *** p<0.01 (the economics convention, same as fixest)
tex_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "^{***}", ifelse(p < 0.05, "^{**}",
         ifelse(p < 0.10, "^{*}", ""))))
}
# Coefficient in math mode so negatives get a real minus sign, e.g. $-1.234^{**}$
tex_coef <- function(b, p, d = 3) sprintf("$%s%s$", formatC(b, format = "f", digits = d),
                                          tex_stars(p))
tex_se   <- function(s, d = 3) sprintf("(%s)", formatC(s, format = "f", digits = d))
tex_n    <- function(n) formatC(n, format = "d", big.mark = ",")

# Writes a booktabs table. head: list of header rows (character vectors).
# rows: list; a character vector becomes "a & b & c \\", a single string starting
# with a backslash (e.g. "\\midrule") is written as-is.
# resize = TRUE shrinks a wide table to the page width (don't use on narrow tables,
# it would blow them up). size is a LaTeX font size such as "\\small".
write_tex_table <- function(file, caption, label, head, rows, notes, align,
                            resize = TRUE, size = NULL) {
  as_line <- function(r) {
    if (length(r) == 1 && startsWith(r, "\\")) r else paste0(paste(r, collapse = " & "), " \\\\")
  }
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    size,
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    if (resize) "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    vapply(head, as_line, ""),
    "\\midrule",
    vapply(rows, as_line, ""),
    "\\bottomrule",
    if (resize) "\\end{tabular}}" else "\\end{tabular}",
    "\\par\\vspace{4pt}",
    sprintf("\\begin{minipage}{\\textwidth}\\footnotesize %s\\end{minipage}", notes),
    "\\end{table}"
  )
  writeLines(lines, file)
}

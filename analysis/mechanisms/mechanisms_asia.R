# Investment mechanism tables, with an "Excluding Asia" sub-panel under each outcome
#
# Produces two tables:
#   1. fdi_panel_asia.tex/csv        -- Anthony's original two treatment types
#                                       (PEPFAR group x Post, HIV prevalence 2003 x Post)
#   2. fdi_panel_triple_asia.tex/csv -- the triple interaction (HIV x PEPFAR x Post)
# Each outcome panel is now two sub-panels: "All countries", then "Excluding Asia",
# same columns throughout -- the same structure the Crown DiD robustness table already
# uses for GDP growth (dropping East, Southeast and South Asia so the 1997-98 Asian
# financial crisis cannot drive the before/after comparison).
#
# East, Southeast and South Asia (World Bank regions East Asia & Pacific and South Asia,
# without the Pacific islands). Central Asia and the Middle East stay in the sample.
# Same list as analysis/replication/twfe_growth_robustness.R, so this table's "Excluding
# Asia" panel is directly comparable to the Crown DiD robustness table's.
# Treated: AFG BGD CHN IDN IND KHM LAO MMR NPL PAK PHL THA TLS VNM (COP: IDN IND KHM MMR VNM)
# Controls: BTN LKA MDV MNG MYS PRK
asia_iso3 <- c("AFG", "BGD", "BTN", "CHN", "IDN", "IND", "KHM", "LAO", "LKA", "MDV",
              "MMR", "MNG", "MYS", "NPL", "PAK", "PHL", "PRK", "THA", "TLS", "VNM")
 
source("analysis/econ_setup.R")
 
out_dir <- "tables/crown_panel/mechanisms"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
 
mech <- read_csv("data/intermediate/wdi_mechanisms_1990_2022.csv", show_col_types = FALSE)
 
mech_outcomes <- c(fdi_pct_gdp     = "FDI, net inflows (\\% of GDP)",
                   fdi_pc          = "FDI, net inflows per person (constant 2022 US\\$)",
                   savings_pct_gdp = "Gross domestic savings (\\% of GDP)",
                   gfcf_pct_gdp    = "Gross fixed capital formation (\\% of GDP)")
transforms <- c(raw = "Raw", win = "Winsorized", winihs = "Winsorized, IHS")
specs      <- c(unadj = "Unadjusted", adj = "Adjusted")
variants   <- c(all = "All countries", excl_asia = "Excluding Asia")
 
panel_m <- left_join(panel, mech, by = c("iso3c", "year"))
for (y in names(mech_outcomes)) {
  w <- winsor(panel_m[[y]])
  panel_m[[paste0(y, "__raw")]]    <- panel_m[[y]]
  panel_m[[paste0(y, "__win")]]    <- w
  panel_m[[paste0(y, "__winihs")]] <- asinh(w)
}
 
# Same make_sample() as everywhere else, with an optional Asia exclusion applied after.
make_sample_variant <- function(group, variant, data = panel_m) {
  d <- make_sample(group, data = data)
  if (variant == "excl_asia") d <- filter(d, !iso3c %in% asia_iso3)
  d
}
 
# =====================================================
# Table 1: Anthony's two treatment types (reuses run_twfe_post from econ_setup.R)
# =====================================================
res1 <- list()
for (g in groups) for (v in names(variants)) {
  d <- make_sample_variant(g, v)
  for (tv in names(treatments)) for (y in names(mech_outcomes))
    for (tr in names(transforms)) for (sp in names(specs)) {
      m  <- run_twfe_post(d, paste0(y, "__", tr), tvar = tv, adjusted = sp == "adj")
      ct <- coeftable(m)[paste0(tv, ":post"), ]
      res1[[length(res1) + 1]] <- tibble(
        outcome = y, group = g, variant = v, treatment = tv, transform = tr, spec = sp,
        est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
        n_obs = nobs(m), n_countries = n_countries(m, d))
    }
}
res1 <- bind_rows(res1)
write_csv(res1, file.path(out_dir, "fdi_panel_asia.csv"))
 
cols1 <- expand.grid(spec = names(specs), group = groups, treatment = names(treatments),
                     stringsAsFactors = FALSE)
pick1 <- function(y, v, tr, j) filter(res1, outcome == y, variant == v, transform == tr,
                                      treatment == cols1$treatment[j],
                                      group == cols1$group[j], spec == cols1$spec[j])
rows1 <- list()
first1 <- TRUE
for (y in names(mech_outcomes)) {
  for (v in names(variants)) {
    if (!first1) rows1[[length(rows1) + 1]] <- "\\midrule"
    first1 <- FALSE
    rows1[[length(rows1) + 1]] <- sprintf(
      "\\multicolumn{9}{l}{\\textit{%s --- %s}} \\\\", mech_outcomes[[y]], variants[[v]])
    for (tr in names(transforms)) {
      r <- lapply(seq_len(nrow(cols1)), \(j) pick1(y, v, tr, j))
      rows1[[length(rows1) + 1]] <- c(transforms[[tr]], sapply(r, \(x) tex_coef(x$est, x$p)))
      rows1[[length(rows1) + 1]] <- c("", sapply(r, \(x) tex_se(x$se)))
    }
    n <- sapply(seq_len(nrow(cols1)), \(j) {
      vv <- filter(res1, outcome == y, variant == v, treatment == cols1$treatment[j],
                   group == cols1$group[j], spec == cols1$spec[j])$n_countries
      if (min(vv) == max(vv)) as.character(vv[1]) else paste0(min(vv), "--", max(vv))
    })
    rows1[[length(rows1) + 1]] <- c("Countries", n)
  }
}
 
write_tex_table(
  file    = file.path(out_dir, "fdi_panel_asia.tex"),
  caption = paste("Investment mechanism: FDI, saving and investment after PEPFAR, two-way",
                  "fixed effects, with and without Asia"),
  label   = "tab:crown_mech_fdi_asia",
  head    = list(c("", sprintf("\\multicolumn{4}{c}{Treatment: %s}", treatments)),
                 "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}",
                 c("", rep(c("\\multicolumn{2}{c}{All PEPFAR sample}",
                             "\\multicolumn{2}{c}{COP sample}"), 2)),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}",
                 c("", rep(c("Unadj.", "Adj."), 4)),
                 c("", sprintf("(%d)", 1:8))),
  rows    = rows1,
  notes   = paste(
    "Each cell is the coefficient on treatment $\\times$ Post from a separate two-way",
    "fixed effects regression with country and year fixed effects; Post = 1 from 2004",
    "for every country. Columns (1)--(4): treatment is the 0/1 PEPFAR (or COP) group",
    "dummy. Columns (5)--(8): treatment is HIV prevalence among ages 15--49 in 2003",
    "(percentage points). Adjusted columns add the 11 baseline (2003) covariates of",
    "Crown et al.\\ (2023), each interacted with Post (columns (6) and (8) leave out",
    "HIV prevalence, which is the treatment there). \\textit{Excluding Asia} drops East,",
    "Southeast and South Asia (20 countries; same list as the Crown DiD robustness",
    "table) so the 1997--98 Asian financial crisis cannot drive the before/after",
    "comparison. Outcomes from the World Bank WDI;",
    "FDI per person is net inflows in current US\\$ deflated to 2022 dollars with the US",
    "GDP deflator and divided by UN WPP population. Winsorized: values outside the",
    "1st--99th percentiles of the full panel set to those percentiles. IHS:",
    "$\\sinh^{-1}(w)$ of the winsorized value. Sample 1990--2022; controls are the 64",
    "comparison countries valid through 2022. Standard errors clustered by country in",
    "parentheses. $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccccccc", resize = TRUE)
 
cat("FDI / investment panel (with Asia robustness) done.\n")
 
# =====================================================
# Table 2: triple interaction (HIV x PEPFAR x Post), Panel A only
# =====================================================
run_twfe_triple <- function(d, y, adjusted = FALSE) {
  rhs <- c("bl_hiv_prev:treat:post", "bl_hiv_prev:post", "treat:post")
  if (adjusted) {
    rhs <- c(rhs, paste0(setdiff(bl_cont, "bl_hiv_prev"), ":post"),
             "i(bl_wb_income, post, ref = 'L')")
  }
  fml <- as.formula(paste(y, "~", paste(rhs, collapse = " + "), "| iso3c + year"))
  feols(fml, data = d, cluster = ~iso3c)
}
 
get_coef <- function(m, vars) {
  ct <- coeftable(m)
  toks <- strsplit(rownames(ct), ":", fixed = TRUE)
  idx <- which(vapply(toks, function(v) setequal(v, vars), logical(1)))
  if (length(idx) != 1) {
    stop("Expected exactly one coefficient matching {", paste(vars, collapse = ", "),
        "} but found ", length(idx), ". Actual coefficient names in this model:\n",
        paste(rownames(ct), collapse = "\n"))
  }
  ct[idx, ]
}
 
res2 <- list()
for (g in groups) for (v in names(variants)) {
  d <- make_sample_variant(g, v)
  for (y in names(mech_outcomes)) for (tr in names(transforms)) for (sp in names(specs)) {
    yy <- paste0(y, "__", tr)
    m  <- run_twfe_triple(d, yy, adjusted = sp == "adj")
    ct <- get_coef(m, c("bl_hiv_prev", "treat", "post"))
    res2[[length(res2) + 1]] <- tibble(
      outcome = y, group = g, variant = v, transform = tr, spec = sp,
      est = ct[["Estimate"]], se = ct[["Std. Error"]], p = ct[["Pr(>|t|)"]],
      n_obs = nobs(m), n_countries = n_countries(m, d))
  }
}
res2 <- bind_rows(res2)
write_csv(res2, file.path(out_dir, "fdi_panel_triple_asia.csv"))
 
cols2 <- expand.grid(spec = names(specs), group = groups, stringsAsFactors = FALSE)
pick2 <- function(y, v, tr, j) filter(res2, outcome == y, variant == v, transform == tr,
                                      group == cols2$group[j], spec == cols2$spec[j])
rows2 <- list()
first2 <- TRUE
for (y in names(mech_outcomes)) {
  for (v in names(variants)) {
    if (!first2) rows2[[length(rows2) + 1]] <- "\\midrule"
    first2 <- FALSE
    rows2[[length(rows2) + 1]] <- sprintf(
      "\\multicolumn{5}{l}{\\textit{%s --- %s}} \\\\", mech_outcomes[[y]], variants[[v]])
    for (tr in names(transforms)) {
      r <- lapply(seq_len(nrow(cols2)), \(j) pick2(y, v, tr, j))
      rows2[[length(rows2) + 1]] <- c(transforms[[tr]], sapply(r, \(x) tex_coef(x$est, x$p)))
      rows2[[length(rows2) + 1]] <- c("", sapply(r, \(x) tex_se(x$se)))
    }
    n <- sapply(seq_len(nrow(cols2)), \(j) {
      vv <- filter(res2, outcome == y, variant == v, group == cols2$group[j],
                   spec == cols2$spec[j])$n_countries
      if (min(vv) == max(vv)) as.character(vv[1]) else paste0(min(vv), "--", max(vv))
    })
    rows2[[length(rows2) + 1]] <- c("Countries", n)
  }
}
 
write_tex_table(
  file    = file.path(out_dir, "fdi_panel_triple_asia.tex"),
  caption = paste("Investment mechanism: HIV prevalence (2003) $\\times$ PEPFAR",
                  "$\\times$ Post, two-way fixed effects, with and without Asia"),
  label   = "tab:crown_mech_fdi_triple_asia",
  head    = list(c("", "\\multicolumn{2}{c}{All PEPFAR sample}",
                       "\\multicolumn{2}{c}{COP sample}"),
                 "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                 c("", rep(c("Unadj.", "Adj."), 2)),
                 c("", sprintf("(%d)", 1:4))),
  rows    = rows2,
  notes   = paste(
    "Each cell is the coefficient on bl\\_hiv\\_prev $\\times$ treat $\\times$ post from",
    "a two-way fixed effects regression with country and year fixed effects; post = 1",
    "from 2004 for every country. Every regression also includes bl\\_hiv\\_prev",
    "$\\times$ post and treat $\\times$ post, required for the triple interaction to be",
    "correctly identified, though only the triple term is reported.",
    "\\textit{Excluding Asia} drops East, Southeast and South Asia (20 countries; same",
    "list as the Crown DiD robustness table) so the 1997--98 Asian financial crisis",
    "cannot drive the before/after comparison.",
    "Adjusted columns add the other 9 continuous baseline (2003) covariates of Crown et",
    "al.\\ (2023), each interacted with post, plus income group $\\times$ post.",
    "Outcomes from the World Bank WDI; FDI per person is net inflows in current US\\$",
    "deflated to 2022 dollars with the US GDP deflator and divided by UN WPP population.",
    "Winsorized: values outside the 1st--99th percentiles of the full panel set to",
    "those percentiles. IHS: $\\sinh^{-1}(w)$ of the winsorized value. Sample",
    "1990--2022; controls are the 64 comparison countries valid through 2022. Standard",
    "errors clustered by country in parentheses.",
    "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  align   = "lcccc", resize = FALSE, size = "\\small")
 
cat("FDI / investment triple-interaction panel (with Asia robustness) done.\n")
cat("Output in", out_dir, "\n")
 
# Crown et al. (2023) 2x2 DiD replication, extended to 2022
# Author: Rishad Vaghaiwalla

library(dplyr)    # tools for editing data
library(readr)    # tools for reading CSV files
library(fixest)   # tools for running regressions

panel <- read_csv("data/intermediate/pepfar_panel_1990_2022_FIXED.csv")

# ---- Build the DiD variables ----
# Time = 1 from 2004 onward (PEPFAR's first year), 0 before
panel <- panel |>
  mutate(Time = as.integer(year >= 2004))
# ---- Regression 1: GDP per capita growth, all PEPFAR, no controls ----
m1 <- feols(gdp_pc_growth ~ Time + grp_pepfar_any + Time:grp_pepfar_any,
            data = panel, cluster = ~iso3c)
summary(m1) 
# ---- Regression 2: GDP per capita growth, COP vs control, no controls ----
# Keep only COP and comparison countries (drop other PEPFAR countries)
panel_cop <- panel |>
  filter(grp_pepfar_cop == 1 | grp_control == 1)

m2 <- feols(gdp_pc_growth ~ Time + grp_pepfar_cop + Time:grp_pepfar_cop,
            data = panel_cop, cluster = ~iso3c)
summary(m2)

# ---- The 11 baseline controls (Crown Table 1), measured before PEPFAR ----
controls <- c("bl_gdp_pc", "bl_pop_total", "bl_life_exp", "bl_fertility",
              "bl_urban_pct", "bl_sec_enroll_gross", "bl_hiv_prev",
              "bl_donor_health_nonpepfar_pc", "bl_dom_hexp_pc_ppp",
              "bl_wb_income", "bl_us_hiv_pre2004")

# ---- Regression 3: GDP per capita growth, all PEPFAR, with controls ----
m3 <- feols(gdp_pc_growth ~ Time + grp_pepfar_any + Time:grp_pepfar_any + .[controls],
            data = panel, cluster = ~iso3c)
summary(m3)
# ---- Regression 4: GDP per capita growth, COP, with controls ----
m4 <- feols(gdp_pc_growth ~ Time + grp_pepfar_cop + Time:grp_pepfar_cop + .[controls],
            data = panel_cop, cluster = ~iso3c)
summary(m4)
# =====================================================
# All outcomes: the same four regressions for each
# =====================================================

# The outcomes from the assignment
outcomes <- c("gdp_pc_growth", "gdp_growth", "gdp",
              "emp_female", "emp_male",
              "prim_disengage_girls", "prim_disengage_boys")

# One treatment column named PEPFAR, so all four models use the same names
panel <- panel |>
  mutate(PEPFAR = grp_pepfar_any)

panel_cop <- panel |>
  filter(grp_pepfar_cop == 1 | grp_control == 1) |>
  mutate(PEPFAR = grp_pepfar_cop)

# A recipe: give it an outcome name (y), get back its four regressions
run_four <- function(y) {
  list(
    "All PEPFAR"     = feols(.[y] ~ Time + PEPFAR + Time:PEPFAR,
                             data = panel, cluster = ~iso3c),
    "All + controls" = feols(.[y] ~ Time + PEPFAR + Time:PEPFAR + .[controls],
                             data = panel, cluster = ~iso3c),
    "COP"            = feols(.[y] ~ Time + PEPFAR + Time:PEPFAR,
                             data = panel_cop, cluster = ~iso3c),
    "COP + controls" = feols(.[y] ~ Time + PEPFAR + Time:PEPFAR + .[controls],
                             data = panel_cop, cluster = ~iso3c)
  )
}
# Run the four regressions for every outcome
results <- list()
for (y in outcomes) {
  results[[y]] <- run_four(y)
}
# Save one results table per outcome (same layout as Crown's tables)
for (y in outcomes) {
  table_text <- capture.output(etable(results[[y]], drop = "bl_"))
  writeLines(table_text, paste0("tables/crown_did_", y, ".txt"))
}


# =====================================================
# LaTeX tables in the Table A8 style (one table per outcome)
# =====================================================
# Each table has four columns, one per regression in run_four():
#   Panel A, All PEPFAR vs. controls: (1) no controls, (2) with the 11 baseline controls
#   Panel B, COP vs. controls:        (3) no controls, (4) with the 11 baseline controls
# Rows: the three DiD coefficients, then group means, then sample sizes.
# Output: tables/crown_did_tex/crown_did_<outcome>.tex, plus crown_did_all.tex,
# which pulls all seven tables into one file for Overleaf.

tex_dir <- "tables/crown_did_tex"
dir.create(tex_dir, recursive = TRUE, showWarnings = FALSE)

# Readable names for the outcomes and coefficients. Percent signs are written
# as \\% because LaTeX treats a bare % as the start of a comment.
outcome_labels <- c(
  gdp_pc_growth        = "GDP per capita growth (\\%)",
  gdp_growth           = "GDP growth (\\%)",
  gdp                  = "GDP (billions of constant 2022 US dollars)",
  emp_female           = "Female employment (\\% of ages 15+)",
  emp_male             = "Male employment (\\% of ages 15+)",
  prim_disengage_girls = "Girls' primary disengagement (share, 0--1)",
  prim_disengage_boys  = "Boys' primary disengagement (share, 0--1)"
)
coef_labels <- c("Time:PEPFAR" = "PEPFAR $\\times$ Post",
                 "PEPFAR"      = "PEPFAR",
                 "Time"        = "Post")

# GDP is stored in dollars, so its coefficients are 12-digit numbers that don't
# print properly in a table. The GDP table uses GDP in billions instead: the same
# regressions, just rescaled.
panel$gdp_bn     <- panel$gdp / 1e9
panel_cop$gdp_bn <- panel_cop$gdp / 1e9

# The data each regression was run on (same order as run_four)
model_data <- list(panel, panel, panel_cop, panel_cop)

# Format a number: commas for big numbers, 3 decimals for small ones
fmt <- function(x) {
  if (is.na(x)) return("--")
  formatC(x, format = "f", digits = if (abs(x) >= 1000) 0 else 3, big.mark = ",")
}

# Group means on the exact rows each regression used (obs(m) = those rows)
#   Control Mean   = comparison countries, 2004 onward
#   Treatment Mean = PEPFAR (or COP) countries, 2004 onward
#   Baseline Mean  = PEPFAR (or COP) countries, before 2004
group_means <- function(m, d, y) {
  est <- d[obs(m), ]
  c(control   = fmt(mean(est[[y]][est$PEPFAR == 0 & est$Time == 1], na.rm = TRUE)),
    treatment = fmt(mean(est[[y]][est$PEPFAR == 1 & est$Time == 1], na.rm = TRUE)),
    baseline  = fmt(mean(est[[y]][est$PEPFAR == 1 & est$Time == 0], na.rm = TRUE)))
}

for (y in outcomes) {
  ms    <- unname(results[[y]])
  y_tab <- y
  if (y == "gdp") {
    ms    <- unname(run_four("gdp_bn"))
    y_tab <- "gdp_bn"
  }
  means <- mapply(group_means, ms, model_data, MoreArgs = list(y = y_tab))
  n_countries <- mapply(function(m, d) length(unique(d$iso3c[obs(m)])), ms, model_data)
  
  etable(ms,
         tex         = TRUE,
         style.tex   = style.tex("aer"),          # clean booktabs look, like Table A8
         digits      = "r3",                      # 3 decimal places
         keep        = c("%^Time:PEPFAR$", "%^PEPFAR$", "%^Time$"),  # hide constant and controls
         order       = c("%^Time:PEPFAR$", "%^PEPFAR$", "%^Time$"),
         dict        = c(coef_labels, outcome_labels),
         headers     = list("^:_:" = list("Panel A: All PEPFAR vs.\\ controls" = 2,
                                          "Panel B: COP vs.\\ controls"        = 2)),
         depvar      = FALSE,
         fitstat     = ~ n,
         signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
         extralines  = list(
           "_Baseline controls" = c("No", "Yes", "No", "Yes"),
           "_Countries"         = n_countries,
           "_Control Mean"      = means["control", ],
           "_Treatment Mean"    = means["treatment", ],
           "_Baseline Mean"     = means["baseline", ]),
         title = paste0("Crown et al.\\ (2023) DiD replication, 1990--2022: ",
                        outcome_labels[[y]]),
         label = paste0("tab:crown_did_", y),
         notes = paste(
           "Each column is a separate pooled OLS regression of the outcome on Post,",
           "PEPFAR and PEPFAR $\\times$ Post (the DiD estimate). Post $= 1$ from 2004.",
           "Panel A compares all PEPFAR countries with all other countries; Panel B",
           "compares COP countries with the Crown comparison group. Columns (2) and (4)",
           "add the 11 baseline covariates of Crown et al.\\ (2023) Table 1 (not shown).",
           "Standard errors clustered by country in parentheses.",
           "Control/Treatment Mean: comparison/treated countries from 2004; Baseline",
           "Mean: treated countries before 2004.",
           "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
         file    = file.path(tex_dir, paste0("crown_did_", y, ".tex")),
         replace = TRUE)
}

# One file that pulls in all seven tables, each on its own page
# (in Overleaf, upload it into the same folder as the seven tables)
writeLines(paste0("\\input{crown_did_", outcomes, "}\n\\clearpage"),
           file.path(tex_dir, "crown_did_all.tex"))

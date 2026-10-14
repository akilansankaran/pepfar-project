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

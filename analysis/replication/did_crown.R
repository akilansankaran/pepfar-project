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
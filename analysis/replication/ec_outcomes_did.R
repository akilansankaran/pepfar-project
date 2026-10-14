suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fixest)
})

panel <- read_csv("data/final/country_panel.csv", show_col_types = FALSE)

controls <- c("gross_cap_form", "trade_open", "inflation", "pop_growth", "primary_enroll")

# The canonical 2x2 form with the three constructed indicators, then the two-way fixed
# effects form in which PEPFAR and Time are absorbed by country and year effects.
m_2x2  <- feols(gdppc_growth ~ PEPFAR + Time + Interaction, panel, cluster = ~iso3)
m_twfe <- feols(gdppc_growth ~ Interaction | iso3 + year, panel, cluster = ~iso3)
m_ctrl <- feols(gdppc_growth ~ Interaction + .[controls] | iso3 + year,
                panel, cluster = ~iso3)
m_focus <- feols(gdppc_growth ~ Interaction_focus + .[controls] | iso3 + year,
                 panel, cluster = ~iso3)

etable(m_2x2, m_twfe, m_ctrl, m_focus,
       file = "tables/replication_did.tex", replace = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fixest)
})

panel <- read_csv("data/final/country_panel.csv", show_col_types = FALSE) |>
  mutate(event_time = ifelse(PEPFAR == 1, pmax(pmin(event_time, 10), -10), -1000))

# Leads are the direct test of pre-trends; the Sun and Abraham estimator guards against
# the contamination of two-way fixed effects under staggered adoption.
es_twfe <- feols(gdppc_growth ~ i(event_time, ref = c(-1, -1000)) | iso3 + year,
                 panel, cluster = ~iso3)
es_sa <- feols(gdppc_growth ~ sunab(ifelse(PEPFAR == 1, pepfar_first_year, 10000), year) |
                 iso3 + year, panel, cluster = ~iso3)

# Country-specific linear trends as the first robustness check on flexible time trends.
m_trend <- feols(gdppc_growth ~ Interaction | iso3[year] + year, panel, cluster = ~iso3)

pdf("tables/event_study.pdf", width = 7, height = 4.5)
iplot(list(es_twfe, es_sa), main = "GDP per capita growth around first PEPFAR disbursement")
dev.off()
etable(es_twfe, es_sa, m_trend, file = "tables/parallel_trends.tex", replace = TRUE)

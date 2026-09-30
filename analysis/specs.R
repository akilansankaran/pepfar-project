# Shared outcome and control sets; analysis scripts should source this rather than
# defining their own lists.

outcomes <- c("gdppc_growth", "gdp_growth", "life_exp", "crude_death_rate")

controls_baseline <- c("gross_cap_form", "trade_open", "inflation", "pop_growth",
                       "primary_enroll")

controls_inflows <- c("fdi_inflow_gdp", "remittances_gdp", "oda_gni", "hiv_aid_non_us")

controls_other_programs <- c()  # HIPC timing and other HIV programs, to be added

controls_full <- c(controls_baseline, controls_inflows, controls_other_programs)

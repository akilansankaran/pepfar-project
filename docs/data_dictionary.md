# Data dictionary

Variables in `data/final/country_panel.csv`. Entries marked TODO are yet to be filled in.

| Variable | Description | Source | Units |
| --- | --- | --- | --- |
| iso3 | ISO 3166 alpha-3 country code | WDI | |
| year | calendar year | | |
| PEPFAR | 1 if the country received any PEPFAR HIV/AIDS disbursement from 2004 | foreignassistance.gov | indicator |
| Time | 1 for years from 2004 | constructed | indicator |
| Interaction | PEPFAR × Time | constructed | indicator |
| PEPFAR_focus | 1 for the fifteen original focus countries | constructed | indicator |
| gdppc_growth | GDP per capita growth | WDI NY.GDP.PCAP.KD.ZG | annual % |
| pepfar_disb_const | US HIV/AIDS disbursements | foreignassistance.gov | constant USD |
| hiv_prev_baseline | HIV prevalence ages 15 to 49 in 2003 | WDI SH.DYN.AIDS.ZS | % |
| TODO | remaining WDI, CRS, GBD, and WPP variables | | |

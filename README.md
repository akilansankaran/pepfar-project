# PEPFAR and Growth Study

This repository holds the data pipeline and analysis code for a study of how the President's Emergency Plan for AIDS Relief (PEPFAR) affected economic growth and individual economic expectations in countries affected by HIV/AIDS. 

The study begins from the country-level evidence of Crown et al. (2023) and Kates et al. (2026), who find growth gains in PEPFAR recipients using
difference-in-differences designs and then proceeds to replicate those estimates on a panel extended with recent
data. We then ask which channels such gains arose from by tying them to the disease burden that the program itself received, and examine the resulting survey evidence on whether individuals' economic expectations shifted once the program began.

Our fundamental research focus is the classification and determination of effects from the introduction of PEPFAR on individual economic expectations and economic growth in
countries affected by HIV/AIDS.

## Components

The **replication** reproduces the published difference-in-differences estimates on a panel of
roughly 157 low and middle income countries, comparing growth in GDP and GDP per capita before and
after 2004 between recipients and non-recipients, across a range of control specifications. The
**parallel trends** module tests the identifying assumption directly through event-study leads and
asks whether the estimates survive flexible country-specific time trends.

The **disease burden** extension interacts the post-program indicator with baseline HIV severity
from the IHME Global Burden of Disease series, so that the growth effect is read against the dose
of disease the program could plausibly relieve. Contemporaneous changes such as HIPC debt relief
and other HIV programs, together with FDI, remittances, and other capital inflows, are controlled
for in order to separate the program from coincident shocks and from reverse causality through
funding intensity.

The **mechanisms** module studies longevity, investment, and individual economic expectations from
the Afrobarometer as channels from the program to growth. **Robustness** is assessed through
data-driven control selection, and the **spillovers** and **counterfactual epidemiology** modules
extend the analysis to regional spillovers and to forward projections of the growth consequences of
ending the epidemic.

## Data

| Source | Content | Access |
| --- | --- | --- |
| World Bank WDI | GDP and GDP per capita growth, controls, FDI, remittances, ODA, HIV prevalence | public API |
| foreignassistance.gov | US HIV/AIDS obligations and disbursements by country and year | public bulk file |
| OECD CRS | HIV/AIDS aid (purpose code 13040) from all donors | public SDMX API |
| IHME GBD | HIV and all-cause deaths, DALYs, prevalence, incidence | public export; annual data restricted |
| UN World Population Prospects | population, life expectancy, under-five mortality, migration | public bulk file |
| World Bank income classification | historical income group by year | public file |
| Afrobarometer | individual assessments of present and expected economic conditions | public, registration required |

## Walkthrough

The repository structure is as follows.
```
data/            raw/ (Dropbox pointer, untracked), interim/, final/
data_pull/       one script per unit of analysis that writes into data/raw
panel_construction/
                 merges the raw pulls into analysis-ready panels in data/final
analysis/        one folder per component, each with a README describing its scope
tables/          regression tables and figures written by the analysis scripts
docs/lit_review/ working bibliography and per-paper notes
setup.R          installs dependencies and records them in renv.lock
```

`data_pull/pull_all_country_data.R` defines one function per country-level source and calls them in
sequence at the bottom, so a single run refreshes every country-year input. GBD exports cannot be
fetched programmatically and are placed by hand in `data/raw/ihme_gbd/`, as described in
`data/raw/README.md`. `data_pull/afrobarometer_pull.R` stays separate because its unit is the
respondent and its coverage is a set of survey rounds over a subset of African countries; it
harmonises the economic-conditions items whose question numbers change across rounds.

`panel_construction/build_panel.R` merges the country pulls, restricts the sample to countries
classified low or middle income at baseline, and constructs the difference-in-differences
variables in one place: `PEPFAR` marks recipients, `Time` marks years from 2004, and `Interaction`
is their product. Alternative treatment definitions (the fifteen original focus countries, funding
per capita, and event time relative to first disbursement) and baseline burden measures are built
alongside them. `panel_construction/build_afrobarometer_panel.R` attaches country treatment status
to respondents and also collapses to country-round means.

## Running

```r
source("setup.R")
source("data_pull/pull_all_country_data.R")
source("data_pull/afrobarometer_pull.R")
source("panel_construction/build_panel.R")
source("panel_construction/build_afrobarometer_panel.R")
source("analysis/replication/did_baseline.R")
source("analysis/parallel_trends/event_study.R")
```

All scripts run from the repository root. Set `PEPFAR_RAW_DIR` to read raw files from the shared
Dropbox folder rather than downloading them again.

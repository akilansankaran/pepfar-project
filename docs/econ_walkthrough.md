# Walkthrough: DiD and event studies for the economic outcomes

Anthony Babu, Oct 2026. This explains every file, every block of code, what data goes in, and what comes out, so you can check it and explain it. Line numbers refer to the files as committed.

---

## 0. The map

Four R files, all run from the repository root (the `pepfar-project` folder):

| Order | File | What it does | Writes to |
|---|---|---|---|
| (shared) | [`analysis/econ_setup.R`](../analysis/econ_setup.R) | Loads the data, defines outcomes, covariates, samples, and helper functions. Not run by itself; the other three start with `source("analysis/econ_setup.R")`, which runs it first. | nothing |
| 1 | [`analysis/replication/did_econ.R`](../analysis/replication/did_econ.R) | **Deliverable 1:** adjusted DiD regressions, All PEPFAR and COP | `tables/crown_did/adjusted/` (`did_adjusted_*.tex`, `did_adjusted.csv`) |
| 2 | [`analysis/parallel_trends/event_study_econ_adjusted.R`](../analysis/parallel_trends/event_study_econ_adjusted.R) | **Deliverable 2:** adjusted event studies, tables and graphs | `tables/crown_event_studies/adjusted/` (`tex/`, `figures/`, `csv/`) |
| 3 | [`analysis/parallel_trends/event_study_economic.R`](../analysis/parallel_trends/event_study_economic.R) | Extra: unadjusted event studies plus robustness (trends, staggered timing, etc.) | `tables/crown_event_studies/robustness/` |

Plus one hand-written LaTeX file, [`tables/crown_results.tex`](../tables/crown_results.tex), which stitches every table and graph into one document. See [`tables/README.md`](../tables/README.md) for the folder map.

To rerun everything (about 20 seconds total):

```
Rscript analysis/replication/did_econ.R
Rscript analysis/parallel_trends/event_study_econ_adjusted.R
Rscript analysis/parallel_trends/event_study_economic.R
```

Packages used: `dplyr`, `readr`, `tidyr` (data handling), `fixest` (regressions), `ggplot2`, `patchwork` (graphs).

---

## 1. The data

**File:** `data/intermediate/pepfar_panel_1990_2022_FIXED.csv`, the team's panel built by Emma/Gemma/Rachel, documented in `data/pepfar_panel_data_dictionary.docx`.

**One row = one country in one year.** 157 countries × 33 years (1990–2022) = 5,181 rows.

Columns this code uses:

| Column | Meaning | How we use it |
|---|---|---|
| `iso3c` | 3-letter country code | identifies countries; SEs are clustered on it; country fixed effects use it |
| `year` | calendar year | year fixed effects, event-study dummies |
| `grp_pepfar_any` | 1 if the country ever got PEPFAR (90 countries) | treatment flag for "All PEPFAR" |
| `grp_pepfar_cop` | 1 if a COP country (31) | treatment flag for "COP" |
| `grp_control_kates` | 1 if a comparison country valid through 2022 (64) | the control group |
| `gdp_pc_growth`, `gdp_growth` | annual % growth (WDI) | outcomes, used as-is |
| `emp_female`, `emp_male` | employment-to-population ratio, % of 15+ (ILO modeled) | outcomes, used as-is |
| `prim_disengage_girls`, `prim_disengage_boys` | out-of-school share of primary-age children (0–1) | outcomes, **×100** to get percentage points like Crown |
| `gdp_pc`, `gdp` | levels in constant 2022 USD | outcomes, **logged** |
| `bl_*` (11 columns) | baseline covariates, fixed at each country's 2003 value | controls in the adjusted models |

**The 11 baseline covariates** (Crown Table 1). Each is one number per country, the same in every year:

| Column | Meaning | Countries missing it (of 154 in the All PEPFAR sample) |
|---|---|---|
| `bl_gdp_pc` | GDP per capita, 2003 | 4 |
| `bl_pop_total` | population | 1 |
| `bl_life_exp` | life expectancy | 2 |
| `bl_fertility` | fertility rate | 2 |
| `bl_urban_pct` | % urban | 2 |
| `bl_sec_enroll_gross` | secondary enrollment | **28** |
| `bl_hiv_prev` | HIV prevalence 15–49 | 3 |
| `bl_donor_health_nonpepfar_pc` | non-PEPFAR donor health aid per person | **24** |
| `bl_dom_hexp_pc_ppp` | domestic health spending per person | 9 |
| `bl_us_hiv_pre2004` | got US HIV aid before 2004 (0/1) | 0 |
| `bl_wb_income` | World Bank income group: L / LM / UM | 6 |

A regression can't use a row with a missing covariate, so **adjusted models only use the 104 countries that have all 11** (98 for employment, 96–101 for disengagement). That's the biggest practical cost of "adjusting for all covariates."

---

## 2. `analysis/econ_setup.R`, block by block

### Lines 10–18: packages
`suppressPackageStartupMessages({ library(...) })` loads the six packages without printing startup chatter. `setFixest_notes(FALSE)` stops fixest from printing "N observations removed because of NA values" for every regression. Those rows are still dropped; it just doesn't announce it.

### Line 20: `ref_year <- 2003`
The reference (omitted) year for event studies: the last year before PEPFAR.

### Lines 26–39: load and prepare the data
```r
panel <- read_csv("data/intermediate/pepfar_panel_1990_2022_FIXED.csv", ...) |>
  mutate(...)
```
`read_csv` reads the file into a table called `panel`. `|>` means "take the thing on the left and feed it into the function on the right." `mutate` adds or changes columns:

- `prim_disengage_girls = 100 * prim_disengage_girls` (and boys): converts 0.23 to 23 so the units are percentage points, matching Crown.
- `log_gdp_pc = log(gdp_pc)`, `log_gdp = log(gdp)`: natural logs. A coefficient of 0.10 on a log outcome means about 10% higher. Without logs, a few huge economies (China, India) would dominate.
- `gdp_bn = gdp / 1e9`: GDP in billions, used only in one robustness check.
- `post = as.integer(year >= 2004)`: Crown's "Time" dummy, 1 from 2004 on, 0 before, **the same for every country**.
- `t = year - ref_year`: a counter that is 0 in 2003, −13 in 1990, +19 in 2022. Used only by the linear-trend models in the robustness script.
- `bl_wb_income = na_if(bl_wb_income, "")`: the CSV stores missing income groups as empty text; this turns them into proper missing values (`NA`) so R drops those rows instead of treating "" as a fourth income group.

### Lines 42–63: outcome names
`outcomes` is a named list: the left side is the column name, the right side is the label for graphs. The order (Crown's five, then our three) controls the column order in every table and graph. `outcomes_tex` holds shorter labels for LaTeX headers, with no `%` signs because LaTeX treats `%` as "start of a comment."

### Lines 67–71: covariate lists
`bl_cont` lists the 10 numeric covariates; `bl_all` adds income group, which is a category and gets handled differently (see below).

### Lines 79–88: samples
```r
make_sample <- function(group, last_year = 2022, control = "grp_control_kates", data = panel) {
  tv <- c("All PEPFAR" = "grp_pepfar_any", "COP" = "grp_pepfar_cop")[[group]]
  data |>
    filter(year <= last_year, .data[[tv]] == 1 | .data[[control]] == 1) |>
    mutate(treat = .data[[tv]])
}
```
Given a group name, this:
1. picks the treatment column (`grp_pepfar_any` or `grp_pepfar_cop`);
2. `filter` keeps rows where the country is in that treatment group **or** is a control (`|` means "or"), up to `last_year`;
3. adds a column `treat` = 1 for treated countries, 0 for controls. Every regression then uses `treat`, so the same code works for both groups.

`samples` runs this for both groups:

| Sample | Rows | Countries | Treated | Controls |
|---|---|---|---|---|
| All PEPFAR | 5,082 | 154 | 90 | 64 |
| COP | 3,135 | 95 | 31 | 64 |

**Dropped entirely:** Venezuela, Colombia, Panama (they got enough PEPFAR money after 2018 to fail the control threshold). In the COP sample, the 59 "other PEPFAR" countries are dropped too: they are neither COP nor controls.

**Decision you could change:** `did_crown.R` (Rishad) keeps Venezuela, Colombia and Panama as controls. I used the Kates-valid group because our window runs to 2022. That's why our numbers differ slightly from his.

`group_file` (line 89) maps "All PEPFAR" → `all_pepfar` and "COP" → `cop` for file names.

### Lines 98–106: `run_es()`, the event-study regression
```r
run_es <- function(data, y, adjusted = FALSE) {
  fe <- if (adjusted) {
    paste0("iso3c + bl_wb_income^year + year[", paste(bl_cont, collapse = ", "), "]")
  } else {
    "iso3c + year"
  }
  f <- as.formula(paste0(y, " ~ i(year, treat, ref = ", ref_year, ") | ", fe))
  feols(f, data = data, cluster = ~iso3c)
}
```
This builds a formula as text, then runs it. For `y = "gdp_pc_growth"` with `adjusted = TRUE`, the formula R actually runs is:

```
gdp_pc_growth ~ i(year, treat, ref = 2003) | iso3c + bl_wb_income^year + year[bl_gdp_pc, bl_pop_total, ..., bl_us_hiv_pre2004]
```

Piece by piece, in fixest syntax:
- **Left of `~`:** the outcome.
- **`i(year, treat, ref = 2003)`:** creates one dummy for each year × `treat`. So `year::1990:treat` is 1 only for treated countries in 1990, and so on. The 2003 dummy is left out. That gives 32 coefficients, `year::1990:treat` … `year::2022:treat`. These are the **β_k** in the graphs.
- **Everything after `|`** is a fixed effect: something fixest absorbs without printing a coefficient.
  - `iso3c`: **country fixed effects**, a separate intercept for each country. This soaks up anything permanent about a country (geography, colonial history, being poor in every year).
  - `bl_wb_income^year`: one fixed effect for **each income group in each year** (3 groups × 33 years = 99). So low-income countries can have their own shock in 2009, lower-middle-income countries a different one, and so on. This also contains ordinary year fixed effects.
  - `year[bl_gdp_pc, ...]`: **varying slopes**. For each of the 10 numeric covariates, it estimates a separate coefficient in every year (10 × 33 = 330). In words: "countries that started with higher HIV prevalence are allowed to follow a different path over time." This is how baseline covariates enter when you have country fixed effects. On their own they'd be absorbed by the country FE, because they never change within a country.
- **`cluster = ~iso3c`:** standard errors clustered by country. A country's errors are allowed to be correlated across its 33 years. This is standard for panels and makes SEs bigger, which is more honest.

Unadjusted (`adjusted = FALSE`) is the same, but the part after `|` is just `iso3c + year`: country FE and year FE.

**What β_k means:** (PEPFAR mean − control mean in year k) − (the same difference in 2003), holding fixed everything the FEs and covariates absorb. If parallel trends holds, β_1990 … β_2002 should all be near zero.

### Lines 109–119: `tidy_es()`
Takes a fitted model `m` and returns a clean table with one row per year: `year`, `est` (coefficient), `se`, `p`, and the 95% CI `lo`/`hi` = est ± 1.96 × se. `coeftable(m)` gives fixest's coefficient table. `grepl("^year::\\d+:treat$", ...)` keeps only the event-study rows; it's a pattern meaning "starts with `year::`, then digits, then `:treat`". `sub(...)` pulls the year number out of the name. 2003 is added back as a row with est = 0 so the graphs show it.

### Lines 121–122: `pre_pat` / `post_pat`
Text patterns that match coefficient names for 1990–2002 and 2004–2022.

### Lines 125–130: `avg_coefs()`
Averages a set of coefficients (e.g. all 19 post-2004 ones) and computes the SE of that average. Math: average = w′b with w = 1/19 for each coefficient; its variance = w′Vw, where V is the coefficients' variance matrix from `vcov(m)`. This accounts for the coefficients being correlated with each other.

### Lines 134–137: `pre_test()`
`wald(m, keep = pre_pat)` runs a **joint Wald (F) test** that all 13 pre-2004 coefficients equal zero at once. It returns a p-value. **Small p (< 0.05) = reject parallel pre-trends.** This is more reliable than eyeballing whether individual CIs cross zero.

### Line 139: `n_countries()`
`obs(m)` lists which rows the regression actually used (after dropping missing data), and this counts the distinct countries among them.

### Lines 141–170: plot style
Colors (`pal`: blue, orange, aqua, from a colorblind-checked palette), text colors, a minimal ggplot theme, and `es_layer()`, which draws the standard event-study look: a horizontal zero line, a dashed vertical line between 2003 and 2004, CI whiskers, and dots. `position_dodge(0.6)` shifts two series slightly left and right so they don't overlap.

### Lines 172–215: LaTeX helpers
- `tex_stars(p)`: `***` if p < 0.01, `**` if < 0.05, `*` if < 0.10 (econ convention, same as fixest).
- `tex_coef(b, p)`: formats a coefficient to 3 decimals inside `$...$` so negatives get a proper minus sign, e.g. `$-3.246^{***}$`.
- `tex_se(s)`: `(0.906)`. `tex_n(n)`: `3,392`.
- `write_tex_table(...)`: assembles a full LaTeX table (`\begin{table}` … `\end{table}`) using `booktabs` rules (`\toprule`, `\midrule`, `\bottomrule`). Each row is a character vector glued together with ` & ` and ended with `\\`. `resize = TRUE` wraps wide tables in `\resizebox{\textwidth}{!}{...}` so they fit the page. `size` sets the font (`\small`). Notes go in a small-print box under the table.

---

## 3. `analysis/replication/did_econ.R`: the DiD (Deliverable 1)

### Lines 26–36: the two regressions

**Panel A, `run_did_crown()`: Crown's own model (pooled OLS).** For `gdp_pc_growth`, R runs exactly:
```
gdp_pc_growth ~ post + treat + post:treat + bl_gdp_pc + bl_pop_total + bl_life_exp + bl_fertility
               + bl_urban_pct + bl_sec_enroll_gross + bl_hiv_prev + bl_donor_health_nonpepfar_pc
               + bl_dom_hexp_pc_ppp + bl_us_hiv_pre2004 + bl_wb_income
```
That is, Y = β0 + β1·Post + β2·PEPFAR + β3·(PEPFAR × Post) + δ′X. No fixed effects. All 33 years are pooled into "before 2004" and "from 2004".
- **β1 (`post`):** how much controls changed after 2004.
- **β2 (`treat`):** how different PEPFAR countries were before 2004.
- **β3 (`post:treat`):** **the DiD estimate**, how much more PEPFAR countries changed than controls.
- Covariates enter as plain levels. `bl_wb_income` is text, so R turns it into dummies (`bl_wb_incomeLM`, `bl_wb_incomeUM`, with L as the base group).

This matches what Rishad's `did_crown.R` does, with all 11 covariates and our control group.

**Panel B, `run_did_twfe()`: two-way fixed effects.** R runs:
```
gdp_pc_growth ~ treat:post | iso3c + bl_wb_income^year + year[bl_gdp_pc, ..., bl_us_hiv_pre2004]
```
This is the adjusted event study from Section 2, but with **one** PEPFAR × Post coefficient instead of 32 year-specific ones. `treat` and `post` don't appear on their own because the country FE and year FE already absorb them. I used the same covariates × year as the event study, so Panel B and the event study are the same model at two levels of detail. (Emma's mortality script uses covariates × Post instead of × year; that's a smaller set of controls.)

### Lines 39–45: `grab()`
Pulls from a model the estimate, SE and p-value for the named terms, plus observations, countries and R². One row per term.

### Lines 47–56: the loop
For each group (2) × outcome (8), it runs both regressions: **16 pooled OLS + 16 TWFE = 32 regressions**. It stacks the results and renames both DiD terms to `"did"`. They're saved to `tables/crown_did/adjusted/did_adjusted.csv`, which is the place to check any number in the tables.

### Lines 62–69: `cell()`
Looks up one number for a table cell, e.g. `cell("COP", "TWFE", "did", "emp_male", "se")` gives `(2.144)`.

### Lines 86–119: separate tables (one per group)
Columns = the 8 outcomes. Panel A rows: PEPFAR × Post, PEPFAR, Post (each with its SE underneath), Observations, Countries, R². Panel B: PEPFAR × Post, Observations, Countries, R². In the COP table the rows read "COP × Post" and "COP".
→ `tables/crown_did/adjusted/did_adjusted_all_pepfar.tex`, `did_adjusted_cop.tex`

### Lines 121–150: combined table
Rows = outcomes, columns = (1) All pooled, (2) All TWFE, (3) COP pooled, (4) COP TWFE. Only the DiD coefficient is shown.
→ `tables/crown_did/adjusted/did_adjusted_combined.tex`

### Results (from `tables/crown_did/adjusted/did_adjusted.csv`)

| Outcome | All: pooled | All: TWFE | COP: pooled | COP: TWFE | Crown (adjusted) |
|---|---|---|---|---|---|
| GDP pc growth | **2.51\*\*\*** | 3.27\* | **2.62\*\*** | 7.27 | 2.07 / 2.50 |
| Female emp. | **−3.25\*\*\*** | 0.28 | **−4.15\*\*\*** | 0.18 | negative, sig. |
| Male emp. | 0.35 | 1.40 | −0.05 | 2.32 | negative, sig. |
| Girls' diseng. | **−10.93\*\*\*** | 1.73 | **−10.75\*\*** | 1.63 | −9.18 / −12.58 |
| Boys' diseng. | **−9.93\*\*\*** | −0.92 | **−12.28\*\*\*** | −5.11 | −7.96 / −12.51 |
| GDP growth | **2.86\*\*\*** | 3.57\* | **2.89\*\*\*** | 7.40 | — |
| log GDP pc | 0.07 | −0.00 | 0.12 | −0.12 | — |
| log GDP | 0.16 | −0.08 | 0.22 | −0.24 | — |

**The pooled OLS column replicates Crown closely** (same signs and similar sizes for GDP growth and disengagement). The remaining differences come from our data running to 2022 and our corrected control group.

**Why does TWFE wipe almost everything out?** In pooled OLS, the covariates only shift the *level* of each country. In TWFE, each covariate gets its own coefficient in every year. So countries with high baseline HIV, low income or low life expectancy are allowed their own path over time. But those characteristics are exactly what made a country a PEPFAR country: baseline HIV prevalence plus pre-2004 US HIV aid alone explain 37% of who is treated. Once each of those gets its own path, the remaining comparison is between PEPFAR and non-PEPFAR countries that looked alike in 2003, and there are few of those. Estimates get much noisier (SEs grow 1.3 to 5 times larger) and move toward zero. **The defensible reading:** the Crown-style result depends on assuming PEPFAR and control countries would have followed the same path despite very different starting points. Once you only compare similar countries, there isn't enough data to say much.

---

## 4. `analysis/parallel_trends/event_study_econ_adjusted.R`: event studies (Deliverable 2)

### Lines 32–50: estimate
For each group × outcome (16 models): `run_es(d, y, adjusted = TRUE)` (Section 2), then
- `tidy_es(m)`: the 32 year coefficients plus the 2003 zero,
- `pre_test(m)`: the joint pre-trend p-value,
- `avg_coefs(m, post_pat)`: the average of the 2004–2022 coefficients and its SE,
- `avg_post_p`: a two-sided p-value for that average from a t distribution with (number of countries − 1) degrees of freedom, the standard choice with clustered SEs.

Saved under `tables/crown_event_studies/adjusted/`: `csv/es_adjusted_coefficients.csv` (every year's coefficient) and `csv/es_adjusted_summary.csv` (one row per model).

### Lines 60–72: `es_cell()` / `s_cell()`
Look up one coefficient/SE, or one summary number, for a table cell. Blank if there's no data that year (e.g. employment starts in 1991).

### Lines 84–123: separate tables (one per group)
Columns = 8 outcomes; rows = years with coefficient and SE underneath. 32 years × 2 lines doesn't fit on one page, so each group gets **two tables**:
- **pre-period (1990–2002)**, ending with the pre-2004 joint test p-value,
- **post-period (2004–2022)**, ending with the average post effect.

Both end with Observations and Countries.
→ `tex/es_adjusted_{all_pepfar,cop}_{pre,post}.tex`

### Lines 125–164: combined tables (one per outcome)
Columns: All PEPFAR coef, SE | COP coef, SE. Rows: every year (2003 marked "reference"), then the pre-test p, the average post effect, observations and countries. Years with no data for either group are dropped.
→ `tex/es_adjusted_combined_<outcome>.tex` (8 files) and `tex/es_adjusted_combined_all.tex`, which inputs all 8, one per page.

### Lines 167–231: graphs
Color always follows the group: **All PEPFAR = blue, COP = orange**. Every graph is saved as **PDF** (sharp, for LaTeX) and **PNG** (for quick viewing). Subtitles print the pre-test p and the average post effect (SE).
- `es_adj_<outcome>_combined`: both groups on one graph (8 graphs)
- `es_adj_<outcome>_all_pepfar`, `es_adj_<outcome>_cop`: each group alone (16 graphs)
- `es_adj_grid_combined`, `es_adj_grid_all_pepfar`, `es_adj_grid_cop`: all 8 outcomes on one page (3 graphs)

How to read one: each dot is β_k with its 95% CI. Left of the dashed line should hug zero if trends were parallel. Right of it is the effect over time.

### Lines 233–270: LaTeX figure blocks
Writes `tex/figures_combined.tex`, `tex/figures_all_pepfar.tex` and `tex/figures_cop.tex`. Each is a list of `\begin{figure} \includegraphics{crown_event_studies/adjusted/figures/....pdf} \caption{...} \end{figure}` blocks. Paths start from the `tables/` folder, and labels start with `crown_` so they don't clash with the mortality tables in the same Overleaf, two per page, so the main document can pull them in with one line each.

### Results (from `csv/es_adjusted_summary.csv`)

| Outcome | All: pre-test p | All: avg post (SE) | COP: pre-test p | COP: avg post (SE) |
|---|---|---|---|---|
| GDP pc growth | 0.022 | −0.21 (2.04) | 0.002 | 1.92 (4.03) |
| Female emp. | 0.016 | 1.14 (1.12) | 0.085 | 1.52 (2.15) |
| Male emp. | 0.312 | 0.34 (1.01) | 0.094 | 0.85 (2.05) |
| Girls' diseng. | 0.015 | −0.85 (2.44) | <0.001 | −3.20 (5.66) |
| Boys' diseng. | 0.133 | −1.99 (2.73) | <0.001 | −5.11 (5.04) |
| GDP growth | 0.046 | 0.19 (2.20) | 0.003 | 2.49 (4.29) |
| log GDP pc | 0.674 | 0.03 (0.09) | <0.001 | 0.16 (0.22) |
| log GDP | 0.624 | −0.01 (0.09) | 0.006 | 0.09 (0.23) |

In words: with all covariates × year, **no average post-2004 effect is significant**, and the pre-trend test still rejects for most outcomes (all of them for COP). The adjusted COP sample is only 53–58 countries with about 400 control parameters, so its early-1990s coefficients swing wildly (e.g. girls' disengagement −38 in 1990). Treat the COP adjusted graphs as very noisy.

---

## 5. `tables/crown_results.tex`
A standard LaTeX article with `booktabs` (table rules), `graphicx` (images) and `amsmath` (equations). It has:
1. a data paragraph,
2. Section 1: the DiD equations, then the combined, All PEPFAR and COP DiD tables,
3. Section 2: the event-study equation, then graphs (combined, All, COP), then tables (combined, All, COP).

Each `\input{...}` pastes in a generated file, by its path from `tables/`. **To compile:** upload `crown_results.tex` with the `crown_did/` and `crown_event_studies/` folders (same folder structure) to Overleaf and compile it. LaTeX isn't installed on this Mac, so I checked the files with a script (balanced braces, matching `\begin`/`\end`, correct column counts, all inputs and figures exist, no stray `%` or `_`), but haven't compiled them. If Overleaf complains, send me the error.

---

## 6. `analysis/parallel_trends/event_study_economic.R` (robustness, from last week)
It now uses the shared setup; its results are unchanged. Sections:
- **2 (lines 35–70):** unadjusted and adjusted event studies for both groups, plus a plain TWFE DiD for comparison (`run_did`).
- **3 (72–112): linear trends.** (a) your slide-3 model, `y ~ treat:post + treat:t | iso3c + year`, where θ = `treat:t` is PEPFAR's extra drift per year and δ = `treat:post` is the jump; (b) country-specific trends, `| iso3c[t] + year`, giving every country its own slope; (c) estimate PEPFAR's 1990–2003 trend, subtract its extension from PEPFAR outcomes, and rerun the event study.
- **4 (114–162): Sun–Abraham.** `sunab(cohort, year)`, where cohort = each country's first PEPFAR year (controls = never). `sa_lead_test()` rebuilds the joint lead test because fixest only reports the averaged coefficients. I checked it reproduces fixest's own numbers exactly.
- **5 (164–206): robustness.** Crown's 1990–2018 window with `grp_control`; drop China and India; countries observed in at least half of the pre and post years; growth winsorized at the 1st/99th percentiles; GDP in levels.
- **6:** CSVs and figures in `tables/crown_event_studies/robustness/`, summarized in its `README.md`.

---

## 7. Decisions worth double-checking with the team
1. **Control group:** `grp_control_kates` (64 countries) rather than Crown's 67. Changing it is one argument in `make_sample()`.
2. **Baseline year:** the panel's `bl_*` covariates are 2003 values; Crown says 2004. 2003 is the cleaner pre-treatment year, but it's a difference from the paper.
3. **Sample window:** 1990–2022 (the team's choice), not Crown's 1990–2018. The robustness script shows the 2018 window gives the same conclusions.
4. **Covariates in TWFE/event study:** × year (flexible) rather than × Post (Emma's mortality script). × year is the honest event-study analog but eats most of the variation.
5. **Missing covariates:** countries missing any of the 11 are dropped (154 → 104). An alternative is fewer covariates or imputing missing ones.
6. **Disengagement units:** percentage points (×100), to match Crown.

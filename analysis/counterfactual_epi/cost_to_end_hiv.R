# Step 6: minimum drug cost of ending HIV over 20 years, lenacapavir only vs. lenacapavir
# plus alimatravir
# Author: Rishad Vaghaiwalla
#
# Assignment notes: price lenacapavir (approved) and alimatravir (Phase 3); treat everyone
# currently living with HIV and keep everyone at risk on PrEP; compare a lenacapavir-only
# world with one where both exist; discount over a 20-year horizon.
#
#   PV cost = sum_{t=0}^{19} (N_treat_t x p_treat + N_prep_t x p_prep_t) / (1 + r)^t
#
# DRUG PRICES ONLY: testing, clinic visits and supply chains come on top, so this is a floor.
# Lenacapavir alone is not a treatment regimen, so treatment is priced at generic TLD.
#
# Two population paths:
#   static   people with HIV and people on PrEP stay at today's numbers for 20 years
#   dynamic  people with HIV evolve as PLHIV_{t+1} = PLHIV_t + new infections_t - AIDS deaths_t,
#            new infections fall linearly to the 2030 target and stay there, and the AIDS death
#            rate among people with HIV stays at its 2025 level. Reports the first year new
#            infections fall below AIDS deaths, Eisinger and Fauci's (2018) criterion for an end.
#
# Every input is in section 0 with its source. Three are ASSUMPTIONS with no source yet
# (alimatravir's price, launch year and share of PrEP users); settle them with the lead.
#
# Run from the repository root:  Rscript analysis/counterfactual_epi/cost_to_end_hiv.R
# Writes to tables/cost_to_end_hiv/{tex, figures, csv}/

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })
source("analysis/econ_setup.R")   # only for the plot theme, palette and LaTeX helpers

# =====================================================
# 0. Inputs
# =====================================================
base_year  <- 2025       # last year of UNAIDS data
start_year <- 2027       # first year of generic lenacapavir supply (expected by end of 2026)
horizon    <- 20         # years, 2027-2046
r_main     <- 0.03       # discount rate, same as the life-years tables
r_alt      <- c(0, 0.05)

# UNAIDS 2026 fact sheet (data for 2025): https://www.unaids.org/resources/fact-sheet
plhiv_0      <- 41.0e6   # people living with HIV
new_inf_0    <- 1.2e6    # new infections
aids_death_0 <- 570e3    # AIDS-related deaths
# UNAIDS (July 2026): about 20 million people need antiretroviral-based prevention
prep_need    <- 20e6
# UN General Assembly Fast-Track target, cited in Eisinger and Fauci (2018): 200,000 new
# infections a year by 2030
target_inf_2030 <- 200e3

# Prices, US$ per person per year (current dollars, no inflation applied)
p_tld <- 45   # generic first-line TLD, Global Fund agreement "below US$45" (Aug 2023)
p_len <- 40   # generic lenacapavir PrEP, Dr. Reddy's / Unitaid / CHAI price for 120 LMICs
              # (FDA tentative approval Nov 2, 2026); US list price is $28,218

# ASSUMPTIONS (no source yet)
p_ali        <- 40    # alimatravir price: placeholder at the lenacapavir price; varied below
ali_launch   <- 2029  # first year of alimatravir supply (Phase 3 readout expected 2027)
ali_share    <- 0.5   # share of PrEP users on alimatravir once it is available
p_ali_grid   <- c(20, 40, 60, 100)   # sensitivity on alimatravir's price

# Optional delivery costs per person-year (0 = drugs only, the minimum the notes ask for)
deliv_treat <- 0
deliv_prep  <- 0

out_dir <- "tables/cost_to_end_hiv"
tex_dir <- file.path(out_dir, "tex"); fig_dir <- file.path(out_dir, "figures"); csv_dir <- file.path(out_dir, "csv")
for (dd in c(tex_dir, fig_dir, csv_dir)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# 1. Population paths, base_year to the end of the horizon
# =====================================================
years <- base_year:(start_year + horizon - 1)
death_rate <- aids_death_0 / plhiv_0

path_static <- tibble(year = years, plhiv = plhiv_0, new_inf = new_inf_0,
                      aids_deaths = aids_death_0, prep_users = prep_need)

new_inf_dyn <- ifelse(years <= 2030,
                      new_inf_0 + (target_inf_2030 - new_inf_0) * (years - base_year) / (2030 - base_year),
                      target_inf_2030)
plhiv_dyn <- numeric(length(years)); plhiv_dyn[1] <- plhiv_0
for (i in seq_along(years)[-1]) {
  plhiv_dyn[i] <- plhiv_dyn[i - 1] + new_inf_dyn[i - 1] - death_rate * plhiv_dyn[i - 1]
}
path_dynamic <- tibble(year = years, plhiv = plhiv_dyn, new_inf = new_inf_dyn,
                       aids_deaths = death_rate * plhiv_dyn, prep_users = prep_need)

paths <- bind_rows(mutate(path_static, path = "Static"), mutate(path_dynamic, path = "Dynamic"))

# =====================================================
# 2. Annual and present-value costs
# =====================================================
annual_costs <- function(paths, sc, p_ali_use = p_ali) {
  with_ali <- sc != "Lenacapavir only"
  paths |>
    filter(year >= start_year) |>
    mutate(scenario  = sc,
           share_ali = if (with_ali) ifelse(year >= ali_launch, ali_share, 0) else 0,
           p_prep    = (1 - share_ali) * p_len + share_ali * p_ali_use,
           cost_treat = plhiv * (p_tld + deliv_treat),
           cost_prep  = prep_users * (p_prep + deliv_prep),
           cost_total = cost_treat + cost_prep,
           t = year - start_year)
}
pv <- function(ac, r) ac |>
  group_by(path, scenario) |>
  summarise(r = r,
            pv_treat = sum(cost_treat / (1 + r)^t), pv_prep = sum(cost_prep / (1 + r)^t),
            pv_total = sum(cost_total / (1 + r)^t), undiscounted_total = sum(cost_total),
            first_year_cost = cost_total[t == 0], last_year_cost = cost_total[t == max(t)],
            .groups = "drop")

ac <- bind_rows(annual_costs(paths, "Lenacapavir only"), annual_costs(paths, "Lenacapavir + alimatravir"))
write_csv(ac, file.path(csv_dir, "annual_costs.csv"))

pv_all <- bind_rows(lapply(c(r_main, r_alt), function(r) pv(ac, r)))
write_csv(pv_all, file.path(csv_dir, "present_value_costs.csv"))

sens <- bind_rows(lapply(p_ali_grid, function(p) annual_costs(paths, "Lenacapavir + alimatravir", p) |>
                           pv(r_main) |> mutate(p_ali = p)))
write_csv(sens, file.path(csv_dir, "sensitivity_alimatravir_price.csv"))

end_year <- path_dynamic |> filter(new_inf < aids_deaths) |> summarise(y = min(year)) |> pull(y)
write_csv(paths, file.path(csv_dir, "population_paths.csv"))

# =====================================================
# 3. Figures
# =====================================================
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
}
p_cost <- ac |>
  filter(path == "Dynamic") |>
  pivot_longer(c(cost_treat, cost_prep), names_to = "part", values_to = "usd") |>
  mutate(part = recode(part, cost_treat = "Treatment (TLD)", cost_prep = "PrEP")) |>
  ggplot(aes(year, usd / 1e9, fill = part)) +
  geom_col(width = 0.8) +
  facet_wrap(~scenario) +
  scale_fill_manual(values = c("Treatment (TLD)" = pal[1], "PrEP" = pal[2])) +
  labs(title = "Annual drug cost of treating everyone with HIV and 20 million people on PrEP",
       subtitle = paste0("Dynamic path; drug prices only. Alimatravir from ", ali_launch, " at $", p_ali,
                         " (placeholder) for ", 100 * ali_share, "% of PrEP users."),
       x = NULL, y = "US$ billion per year") +
  theme_es
save_both(p_cost, "annual_drug_cost", 9, 4.5)

p_epi <- path_dynamic |>
  select(year, `New infections` = new_inf, `AIDS deaths` = aids_deaths) |>
  pivot_longer(-year) |>
  ggplot(aes(year, value / 1e6, color = name)) +
  geom_line(linewidth = 0.9) +
  { if (is.finite(end_year)) geom_vline(xintercept = end_year, color = ink2, linetype = "dashed", linewidth = 0.4) } +
  scale_color_manual(values = c("New infections" = pal[2], "AIDS deaths" = pal[1])) +
  labs(title = "Dynamic path: new infections vs. AIDS deaths",
       subtitle = if (is.finite(end_year)) paste0("Infections fall below deaths in ", end_year,
                                                 " (dashed), the Eisinger-Fauci criterion for an end") else
                   "Infections stay above deaths over the horizon",
       x = NULL, y = "Millions per year") +
  theme_es
save_both(p_epi, "dynamic_infections_vs_deaths", 8, 4)

# =====================================================
# 4. LaTeX table
# =====================================================
fmt_b <- function(v) sprintf("%.1f", v / 1e9)
row_for <- function(pth, sc) {
  q <- function(r) filter(pv_all, path == pth, scenario == sc, r == !!r)
  c(paste0(pth, ": ", sc), fmt_b(q(r_main)$pv_treat), fmt_b(q(r_main)$pv_prep), fmt_b(q(r_main)$pv_total),
    fmt_b(q(0)$pv_total), fmt_b(q(0.05)$pv_total), fmt_b(q(r_main)$first_year_cost))
}
rows <- list(row_for("Static", "Lenacapavir only"), row_for("Static", "Lenacapavir + alimatravir"),
             row_for("Dynamic", "Lenacapavir only"), row_for("Dynamic", "Lenacapavir + alimatravir"))
write_tex_table(
  file    = file.path(tex_dir, "cost_to_end_hiv.tex"),
  caption = paste0("Minimum drug cost of treating everyone with HIV and keeping those at risk on PrEP, ",
                   start_year, "--", start_year + horizon - 1),
  label   = "tab:cost_to_end_hiv",
  head    = list(c("", "\\multicolumn{3}{c}{Present value, 3\\%}", "PV, 0\\%", "PV, 5\\%", "First-year"),
                 "\\cmidrule(lr){2-4}",
                 c("Path: scenario", "Treatment", "PrEP", "Total", "Total", "Total", "cost")),
  rows    = rows,
  notes   = paste(
    "US\\$ billion, drug prices only (current dollars, no inflation). Treatment: everyone living with HIV on",
    paste0("generic TLD at \\$", p_tld, " per person-year. PrEP: ", prep_need / 1e6, " million people (UNAIDS) on generic"),
    paste0("lenacapavir at \\$", p_len, "; in the second scenario ", 100 * ali_share, "\\% switch to alimatravir from ",
           ali_launch, " at a placeholder \\$", p_ali, ","),
    "with no price set yet. Static: 2025 numbers held constant (41.0 million people with HIV). Dynamic:",
    "new infections fall linearly from 1.2 million (2025) to the 200,000 Fast-Track target in 2030, and the",
    paste0("2025 AIDS death rate among people with HIV is held constant; infections fall below deaths in ",
           ifelse(is.finite(end_year), end_year, "no year of the horizon"), "."),
    "Delivery costs (testing, visits, supply chain) are excluded, so these are lower bounds."),
  align = "lcccccc", resize = FALSE)

# =====================================================
# 5. Console summary
# =====================================================
cat("\nPresent value of drug costs, US$ billion:\n")
print(pv_all |> mutate(across(c(pv_treat, pv_prep, pv_total, undiscounted_total, first_year_cost, last_year_cost),
                              ~ round(.x / 1e9, 1))), n = Inf, width = Inf)
cat("\nAlimatravir price sensitivity (dynamic and static, 3%):\n")
print(sens |> transmute(path, p_ali, pv_total_bn = round(pv_total / 1e9, 1)), n = Inf)
cat("\nDynamic path: new infections fall below AIDS deaths in", end_year, "\n")
cat("\nDone. Output in", out_dir, "\n")

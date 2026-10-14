# Kates panel replication: country + year FE, shared Post (year >= 2004)
#   Y_it = a_i + g_t + b (PEPFAR_i x Post_t) [+ X_i^bl x Post_t] + e_it
# 3 outcomes x 3 specs x 2 samples (PEPFAR-all / COP) = 18 regressions

# install.packages(c("fixest", "modelsummary", "dplyr"))
library(fixest)
library(modelsummary)
library(dplyr)

d <- read.csv("~/Desktop/pepfar_panel_1990_2022_FIXED.csv") |>
  filter(year >= 1990, year <= 2022) |>
  mutate(
    post         = as.integer(year >= 2004),
    bl_wb_income = na_if(bl_wb_income, "")
  )

# Baseline covariates (time-invariant; absorbed by country FE on their own,
# so they enter interacted with Post)
bl_1_9  <- c("bl_gdp_pc", "bl_pop_total", "bl_life_exp", "bl_fertility",
             "bl_urban_pct", "bl_sec_enroll_gross", "bl_hiv_prev",
             "bl_donor_health_nonpepfar_pc", "bl_dom_hexp_pc_ppp")
bl_1_11 <- c(bl_1_9, "bl_wb_income", "bl_us_hiv_pre2004")

outcomes <- c(mort_allcause = "All-cause mortality",
              mort_maternal = "Maternal mortality",
              mort_child    = "Child mortality")

# Samples: treated group vs. Kates-valid controls (grp_control_kates drops the
# 3 controls that exceed the $1M / $0.05pc PEPFAR-linked disbursement threshold)
samples <- list(
  "PEPFAR-all" = d |> filter(grp_pepfar_any == 1 | grp_control_kates == 1) |>
                      mutate(treat = grp_pepfar_any),
  "COP"        = d |> filter(grp_pepfar_cop == 1 | grp_control_kates == 1) |>
                      mutate(treat = grp_pepfar_cop)
)

make_fml <- function(y, ctrls) {
  rhs <- "treat:post"
  if (length(ctrls)) {
    # factor baseline (income group) via i(); numeric ones as x:post
    num <- setdiff(ctrls, "bl_wb_income")
    rhs <- c(rhs, paste0(num, ":post"))
    if ("bl_wb_income" %in% ctrls) rhs <- c(rhs, "i(bl_wb_income, post, ref = 'L')")
  }
  as.formula(paste(y, "~", paste(rhs, collapse = " + "), "| iso3c + year"))
}

specs <- list("Unadj." = character(0), "BL 1-9" = bl_1_9, "BL 1-11" = bl_1_11)

models <- list()
for (s in names(samples)) for (y in names(outcomes)) for (sp in names(specs)) {
  models[[s]][[y]][[sp]] <- feols(make_fml(y, specs[[sp]]),
                                  data = samples[[s]], cluster = ~iso3c)
}

# ---- Tables: one per sample, columns = outcome x spec ----------------------
gof_map <- list(
  list(raw = "nobs",      clean = "Observations", fmt = 0),
  list(raw = "r.squared", clean = "R²",           fmt = 3)
)

n_countries <- function(ms) sapply(ms, \(m) length(unique(m$fixef_id$iso3c)))

make_table <- function(s) {
  ms <- unlist(models[[s]], recursive = FALSE)
  names(ms) <- rep(names(specs), length(outcomes))
  rows <- data.frame(
    term = c("Country FE", "Year FE", "Baseline controls × Post", "Countries"),
    rbind(rep("Yes", length(ms)), rep("Yes", length(ms)),
          rep(c("None", "Vars 1–9", "Vars 1–11"), length(outcomes)),
          n_countries(ms))
  )
  attr(rows, "position") <- 3:6
  modelsummary(
    ms,
    output    = "tinytable",
    coef_map  = c("treat:post" = "PEPFAR × Post"),
    stars     = c("*" = .1, "**" = .05, "***" = .01),
    gof_map   = gof_map,
    add_rows  = rows,
    title     = paste0("PEPFAR × Post (TWFE), ", s, " vs. controls, 1990–2022"),
    notes     = "SEs clustered by country. Post = 1 if year ≥ 2004 for all countries."
  ) |>
    tinytable::group_tt(j = setNames(list(2:4, 5:7, 8:10), outcomes))
}

tab_all <- make_table("PEPFAR-all")
tab_cop <- make_table("COP")
tab_all
tab_cop

# Save (format from extension: .html, .docx, .tex, .png, ...)
tinytable::save_tt(tab_all, "kates_panel_pepfar_all.html", overwrite = TRUE)
tinytable::save_tt(tab_cop, "kates_panel_cop.html",        overwrite = TRUE)

# parallel_trends

Tests the identifying assumption behind the replication. `event_study.R` estimates leads and lags
around each country's first PEPFAR disbursement, both by two-way fixed effects and by the Sun and
Abraham interaction-weighted estimator, and refits the headline model with country-specific trends.

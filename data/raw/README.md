# data/raw

Raw downloads are not tracked in git. They live in the shared project Dropbox folder, and the
scripts in `data_pull/` write into this directory by default. Set the environment variable
`PEPFAR_RAW_DIR` to point the pull scripts at the Dropbox copy instead of re-downloading.

Sources that cannot be fetched programmatically (the IHME GBD Results Tool exports and the
Afrobarometer merged-round files, which require registration) should be placed here by hand under
the file names expected in `data_pull/pull_all_country_data.R` and `data_pull/afrobarometer_pull.R`.
Restricted annual IHME data, if granted, must stay in Dropbox and never enter the repository.

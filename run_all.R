# Runs the R pipeline (R/01-05). Expects the A4/LEARN public-release CSVs in data/raw/.
# Usage from the repo root:  Rscript run_all.R > output/run_all.log 2>&1
# The SAS step and the reconciliation (R/06, sas/, R/07) are run separately; see README.
for (s in sprintf("R/%02d_%s.R", 1:5, c("derive_adsl", "derive_adqs", "mmrm_primary", "amyloid_slope", "tables_figures")))
  { message("== ", s); source(s, echo = FALSE) }
cat("\n"); print(sessionInfo())

# Runs the full pipeline. Expects the A4/LEARN public-release CSVs in data/raw/.
for (s in sprintf("R/%02d_%s.R", 1:5, c("derive_adsl", "derive_adqs", "mmrm_primary", "amyloid_slope", "tables_figures")))
  { message("== ", s); source(s, echo = FALSE) }

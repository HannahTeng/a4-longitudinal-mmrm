# Reconcile R (mmrm) and SAS (PROC MIXED) treatment-difference estimates after running sas/01_mmrm_primary.sas
source("R/00_config.R")
r   <- read_csv(file.path(OUT_TAB, "aim1_mmrm_tx_difference.csv"), show_col_types = FALSE) %>%
  transmute(AVISIT, r_est = estimate, r_se = SE, r_df = df)
sas <- read_csv(file.path(OUT_TAB, "aim1_sas_tx_difference.csv"), show_col_types = FALSE) %>%
  transmute(AVISIT, sas_est = Estimate, sas_se = StdErr, sas_df = DF)
rec <- inner_join(r, sas, by = "AVISIT") %>%
  mutate(d_est = abs(r_est - sas_est), d_se = abs(r_se - sas_se), d_df = abs(r_df - sas_df))
print(rec); write_csv(rec, file.path(OUT_TAB, "aim1_reconciliation_r_vs_sas.csv"))
cat(sprintf("Max |diff| estimate = %.2e, SE = %.2e, df = %.2f\n", max(rec$d_est), max(rec$d_se), max(rec$d_df)))

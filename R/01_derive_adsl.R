# ------------------------------------------------------------------
# 01_derive_adsl.R — subject-level analysis dataset (ADSL) from SUBJINFO
# Input : data/raw/SUBJINFO.csv  (A4/LEARN public release, Synapse syn61250768)
# Output: data/derived/adsl.rds, output/tables/adsl_summary.csv
# ------------------------------------------------------------------
source("R/00_config.R")

subj <- read_raw("SUBJINFO")

adsl <- subj %>%
  transmute(
    BID,
    SUBSTUDY  = factor(SUBSTUDY, levels = c("A4", "LEARN", "SF")),
    TX        = factor(TX, levels = c("Placebo", "Solanezumab")),
    # --- flags (as released) ---
    MITTFL, SUBJITTTR, SUBJSAFTR, LEARNFL, TAUPETFL, PCCPFL, SUBJCOMPTR,
    # --- baseline covariates for the MMRM ---
    AGEYR     = as.numeric(AGEYR),
    EDCCNTU   = as.numeric(EDCCNTU),
    SEX       = factor(SEX, levels = c(1, 2), labels = c("Female", "Male")),   # per derived_datadic
    APOE4     = factor(APOEGNPRSNFLG, levels = c(0, 1), labels = c("non-carrier", "carrier")),
    APOE4_AN  = factor(AAPOEGNPRSNFLG, levels = c(0, 1), labels = c("non-carrier", "carrier")),  # analysis flag: missing -> non-carrier (randomization plan)
    SUVRCER   = as.numeric(SUVRCER),         # florbetapir SUVR (cerebellar reference)
    AMYLCENT  = as.numeric(AMYLCENT),        # centiloids
    MMSEBL    = as.numeric(MMSETSV1),
    PACC_W240_SUBJ = as.numeric(PACCV6),     # subject-level PACC at V6; cross-checked vs ADQS later
    TRTDUR_D  = as.numeric(TRTENDT_DAYS_T0) - as.numeric(TRTSTDT_DAYS_T0),
    FRAMINGHAM
  ) %>%
  mutate(
    AMYL_TERT = if_else(SUBSTUDY %in% c("A4", "LEARN"),
                        ntile(AMYLCENT, 3) |> factor(labels = c("T1 (low)", "T2", "T3 (high)")),
                        factor(NA)),
    COHORT    = case_when(SUBSTUDY == "A4"    ~ paste("A4", TX),
                          SUBSTUDY == "LEARN" ~ "LEARN (Abeta-)",
                          TRUE ~ "Screen failure") |> factor()
  )

# --- derivation QC: every randomized subject must have a TX; mITT subset of randomized
stopifnot(all(!is.na(adsl$TX[adsl$SUBSTUDY == "A4"])))
stopifnot(all(adsl$SUBJITTTR[adsl$MITTFL == 1] == 1))

saveRDS(adsl, file.path(DERIVED_DIR, "adsl.rds"))

# --- Table 1-style baseline summary (by cohort), written as CSV for the report
adsl %>%
  filter(SUBSTUDY != "SF") %>%
  group_by(COHORT) %>%
  summarise(
    n            = n(),
    age_mean_sd  = sprintf("%.1f (%.1f)", mean(AGEYR, na.rm = TRUE), sd(AGEYR, na.rm = TRUE)),
    edu_mean_sd  = sprintf("%.1f (%.1f)", mean(EDCCNTU, na.rm = TRUE), sd(EDCCNTU, na.rm = TRUE)),
    apoe4_pct    = sprintf("%.1f%%", 100 * mean(APOE4 == "carrier", na.rm = TRUE)),
    centiloid    = sprintf("%.1f (%.1f)", mean(AMYLCENT, na.rm = TRUE), sd(AMYLCENT, na.rm = TRUE)),
    mmse_bl      = sprintf("%.1f (%.1f)", mean(MMSEBL, na.rm = TRUE), sd(MMSEBL, na.rm = TRUE)),
    completed_pct= sprintf("%.1f%%", 100 * mean(SUBJCOMPTR == 1, na.rm = TRUE))
  ) %>%
  write_csv(file.path(OUT_TAB, "adsl_summary.csv"))

message("ADSL: ", nrow(adsl), " subjects; mITT n = ", sum(adsl$MITTFL == 1),
        "; OBS (placebo + LEARN) n = ", nrow(POP_OBS(adsl)))

# ------------------------------------------------------------------
# 02_derive_adqs.R — longitudinal analysis dataset (ADQS-style) for PACC
# Input : data/raw/PACC.csv, data/derived/adsl.rds
# Output: data/derived/adqs_pacc.rds
# One row per subject x scheduled double-blind visit (baseline + every 24 wk to wk 240).
# ------------------------------------------------------------------
source("R/00_config.R")

adsl <- readRDS(file.path(DERIVED_DIR, "adsl.rds"))
pacc <- read_raw("PACC") %>% rename(PACC_RAW = `PACC.raw`, PACC_CHG = PACC) %>%
  mutate(VISCODE = suppressWarnings(as.integer(VISCODE)))

# Scheduled visits carrying the PACC: VISCODE 6 (baseline), 12, 18, ..., 66 (wk 240)
# Week = (VISCODE - 6) * 4.  Codes > 66 are open-label; 701/997/999 are unscheduled/ET.
adqs <- pacc %>%
  filter(SUBSTUDY %in% c("A4", "LEARN"), VISCODE >= 6, VISCODE <= 66, VISCODE %% 6 == 0) %>%
  mutate(
    AVISITN = (VISCODE - 6) * 4,
    AVISIT  = factor(paste0("Week ", AVISITN), levels = paste0("Week ", PACC_VISIT_WEEKS)),
    ABLFL   = if_else(VISCODE == 6, "Y", "")
  ) %>%
  inner_join(adsl %>% select(BID, TX, COHORT, MITTFL, AGEYR, EDCCNTU, SEX, APOE4, APOE4_AN, SUVRCER, AMYLCENT, MMSEBL),
             by = "BID") %>%
  group_by(BID) %>%
  mutate(BASE = PACC_RAW[AVISITN == 0][1],
         CHG  = PACC_RAW - BASE) %>%            # recomputed; cross-checked below vs released PACC_CHG
  ungroup() %>%
  select(BID, SUBSTUDY, COHORT, TX, MITTFL, AVISIT, AVISITN, ABLFL,
         AVAL = PACC_RAW, BASE, CHG, PACC_CHG_RELEASED = PACC_CHG, PACCIMPFLG,
         AGEYR, EDCCNTU, SEX, APOE4, APOE4_AN, SUVRCER, AMYLCENT, MMSEBL) %>%
  arrange(BID, AVISITN)

# --- derivation QC
chk <- adqs %>% filter(!is.na(CHG), !is.na(PACC_CHG_RELEASED)) %>%
  summarise(max_abs_diff = max(abs(CHG - PACC_CHG_RELEASED)))
message(sprintf("CHG vs released change-from-baseline: max |diff| = %.2e", chk$max_abs_diff))
stopifnot(chk$max_abs_diff < 1e-5)
stopifnot(!any(duplicated(adqs[, c("BID", "AVISITN")])))

saveRDS(adqs, file.path(DERIVED_DIR, "adqs_pacc.rds"))

adqs %>% count(COHORT, AVISIT) %>% pivot_wider(names_from = COHORT, values_from = n) %>%
  write_csv(file.path(OUT_TAB, "pacc_visit_completion.csv"))
message("ADQS rows: ", nrow(adqs), " | subjects: ", n_distinct(adqs$BID))

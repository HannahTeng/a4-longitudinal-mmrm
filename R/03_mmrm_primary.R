# ------------------------------------------------------------------
# 03_mmrm_primary.R — Aim 1: A4 primary analysis replication
# MMRM on change from baseline in PACC (mITT, double-blind period, wk 24-240)
#   CHG ~ TX*AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN, unstructured covariance, Kenward-Roger df
# Compared against the published A4 result (Sperling et al., NEJM 2023).
# ------------------------------------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({ library(mmrm); library(emmeans); library(nlme) })

adqs <- readRDS(file.path(DERIVED_DIR, "adqs_pacc.rds"))

mitt <- adqs %>%
  filter(SUBSTUDY == "A4", MITTFL == 1, AVISITN > 0, !is.na(CHG), !is.na(BASE)) %>%
  filter(!is.na(AGEYR), !is.na(EDCCNTU), !is.na(APOE4_AN)) %>%
  mutate(AVISIT = droplevels(AVISIT), BID = factor(BID), TX = relevel(TX, "Placebo"))

message("MMRM analysis set: ", n_distinct(mitt$BID), " subjects, ", nrow(mitt), " post-baseline PACC records")

fit_un <- mmrm(CHG ~ TX * AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN + us(AVISIT | BID),
               data = mitt, reml = TRUE, method = "Kenward-Roger")

# LS-means by arm and visit; treatment difference at each visit (primary: Week 240)
em   <- emmeans(fit_un, ~ TX | AVISIT)
lsm  <- as.data.frame(em)
diff <- as.data.frame(contrast(em, method = "revpairwise")) %>%   # Solanezumab - Placebo
  mutate(lower = estimate - qt(0.975, df) * SE, upper = estimate + qt(0.975, df) * SE)

write_csv(lsm,  file.path(OUT_TAB, "aim1_mmrm_lsmeans.csv"))
write_csv(diff, file.path(OUT_TAB, "aim1_mmrm_tx_difference.csv"))
saveRDS(fit_un, "output/aim1_mmrm_fit.rds")

cat("\n=== Treatment difference (Solanezumab - Placebo), change from baseline PACC ===\n")
print(diff %>% select(AVISIT, estimate, SE, df, lower, upper, p.value) %>%
        mutate(across(c(estimate, SE, lower, upper), ~ round(.x, 3)), p.value = signif(p.value, 2)))

# --- QC cross-check: same model via nlme::gls (independent implementation)
gls_fit <- gls(CHG ~ TX * AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN,
               data = mitt, correlation = corSymm(form = ~ as.integer(AVISIT) | BID),
               weights = varIdent(form = ~ 1 | AVISIT), method = "REML",
               control = glsControl(maxIter = 200, msMaxIter = 200))
cf_mm  <- coef(fit_un); cf_gls <- coef(gls_fit)
qc <- tibble(term = names(cf_mm), mmrm = cf_mm, gls = cf_gls[names(cf_mm)]) %>%
  mutate(abs_diff = abs(mmrm - gls))
write_csv(qc, file.path(OUT_TAB, "aim1_qc_mmrm_vs_gls.csv"))
cat(sprintf("\nQC: max |coef(mmrm) - coef(gls)| = %.2e\n", max(qc$abs_diff)))

# ------------------------------------------------------------------
# 04_amyloid_slope.R — Aim 2: amyloid burden x time on PACC in untreated participants
# Population: A4 placebo arm + LEARN (amyloid-negative), double-blind-period visits
# (a) random-intercept/slope LMM:  PACC ~ AMYLCENT*WEEK + covariates + (1 + WEEK | BID)
# (b) MMRM with amyloid tertile x visit (categorical-time analogue)
# (c) covariance-structure comparison for the marginal model (UN / AR1 / CS / heterogeneous)
# ------------------------------------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({ library(lme4); library(nlme); library(mmrm); library(emmeans) })

adqs <- readRDS(file.path(DERIVED_DIR, "adqs_pacc.rds"))
obs <- adqs %>%
  filter(COHORT %in% c("A4 Placebo", "LEARN (Abeta-)"), !is.na(AVAL), !is.na(AMYLCENT),
         !is.na(AGEYR), !is.na(EDCCNTU), !is.na(APOE4_AN)) %>%
  mutate(YEARS = AVISITN / 52.18, CENT10 = AMYLCENT / 10,          # centiloids per 10 units
         BID = factor(BID), AVISIT = droplevels(AVISIT))
message("Aim 2 analysis set: ", n_distinct(obs$BID), " subjects, ", nrow(obs), " PACC records")

# (a) random slopes
lmm <- lmer(AVAL ~ CENT10 * YEARS + AGEYR + EDCCNTU + APOE4_AN + (1 + YEARS | BID), data = obs, REML = TRUE)
cs <- summary(lmm)$coefficients
lmm_tab <- tibble(term = rownames(cs), estimate = cs[, "Estimate"], std.error = cs[, "Std. Error"], statistic = cs[, "t value"])
write_csv(lmm_tab, file.path(OUT_TAB, "aim2_lmm_random_slope.csv"))
cat("\n=== Random-slope LMM fixed effects (PACC units) ===\n"); print(lmm_tab, n = 20)
cat(sprintf("Interpretation: each +10 centiloids changes the annual PACC slope by %.3f (SE %.3f)\n",
            lmm_tab$estimate[lmm_tab$term == "CENT10:YEARS"], lmm_tab$std.error[lmm_tab$term == "CENT10:YEARS"]))

# (b) MMRM by amyloid tertile (post-baseline change)
obs_chg <- obs %>% filter(AVISITN > 0, !is.na(CHG)) %>%
  mutate(AMYL_GRP = cut(AMYLCENT, c(-Inf, 20, 50, Inf), labels = c("<20 CL", "20-50 CL", ">50 CL")))
fit_grp <- mmrm(CHG ~ AMYL_GRP * AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN + us(AVISIT | BID),
                data = obs_chg, method = "Kenward-Roger")
lsm_grp <- as.data.frame(emmeans(fit_grp, ~ AMYL_GRP | AVISIT))
write_csv(lsm_grp, file.path(OUT_TAB, "aim2_mmrm_lsmeans_by_amyloid.csv"))

# (c) covariance-structure comparison (gls, same mean model, REML -> compare AIC/BIC)
mm <- CHG ~ AMYL_GRP * AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN
tvis <- ~ as.integer(AVISIT) | BID
cands <- list(
  "Compound symmetry"            = gls(mm, obs_chg, correlation = corCompSymm(form = tvis), method = "REML"),
  "AR(1)"                        = gls(mm, obs_chg, correlation = corAR1(form = tvis),      method = "REML"),
  "AR(1), heterogeneous var"     = gls(mm, obs_chg, correlation = corAR1(form = tvis),
                                       weights = varIdent(form = ~ 1 | AVISIT), method = "REML"),
  "Unstructured"                 = gls(mm, obs_chg, correlation = corSymm(form = tvis),
                                       weights = varIdent(form = ~ 1 | AVISIT), method = "REML",
                                       control = glsControl(maxIter = 200, msMaxIter = 200))
)
cov_tab <- tibble(structure = names(cands),
                  n_cov_params = sapply(cands, function(f) length(unlist(f$modelStruct))),
                  AIC = sapply(cands, AIC), BIC = sapply(cands, BIC), logLik = sapply(cands, logLik)) %>%
  arrange(AIC)
lrt_p <- anova(cands[["AR(1), heterogeneous var"]], cands[["Unstructured"]])$`p-value`[2]
write_csv(cov_tab, file.path(OUT_TAB, "aim2_covariance_comparison.csv"))
cat("\n=== Covariance structure comparison (REML) ===\n"); print(cov_tab)
cat(sprintf("LRT heterogeneous AR(1) vs unstructured: p = %.3g\n", lrt_p))
saveRDS(list(lmm = lmm, mmrm_grp = fit_grp), "output/aim2_fits.rds")

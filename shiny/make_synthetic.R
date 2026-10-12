# ------------------------------------------------------------------
# make_synthetic.R: simulate a longitudinal cognitive-study dataset with the
# STRUCTURE of A4/LEARN (3 cohorts, PACC every 24 weeks to week 240, MAR-type dropout)
# from the fitted Aim 2 random-slope model, so the public dashboard contains no
# real participant data.  Run from the repo root after R/04_amyloid_slope.R.
# ------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
set.seed(20261012)
fits <- readRDS("output/aim2_fits.rds"); lmm <- fits$lmm
beta <- lme4::fixef(lmm); vc <- as.data.frame(lme4::VarCorr(lmm))
sd_int <- vc$sdcor[vc$grp == "BID" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
sd_slp <- vc$sdcor[vc$grp == "BID" & vc$var1 == "YEARS" & is.na(vc$var2)]
rho    <- vc$sdcor[vc$grp == "BID" & !is.na(vc$var2)]
sd_res <- vc$sdcor[vc$grp == "Residual"]
weeks  <- seq(0, 240, by = 24)

n <- c("Placebo" = 220, "Active" = 220, "Amyloid-negative" = 160)
subj <- bind_rows(lapply(names(n), function(g) tibble(
  COHORT = g, AGEYR = round(rnorm(n[g], if (g == "Amyloid-negative") 70.5 else 72, 4.7), 1),
  EDCCNTU = pmin(20, pmax(8, round(rnorm(n[g], 16.6, 2.8)))),
  SEX = sample(c("Female", "Male"), n[g], TRUE, c(0.6, 0.4)),
  APOE4 = sample(c("non-carrier", "carrier"), n[g], TRUE, if (g == "Amyloid-negative") c(0.77, 0.23) else c(0.41, 0.59)),
  AMYLCENT = if (g == "Amyloid-negative") round(rnorm(n[g], 4, 9), 1) else round(pmax(20, rnorm(n[g], 66, 32)), 1)
))) %>% mutate(BID = sprintf("S%04d", row_number()), SITE = sample(sprintf("Site %02d", 1:12), n(), TRUE))

# random effects
re <- MASS::mvrnorm(nrow(subj), c(0, 0), matrix(c(sd_int^2, rho*sd_int*sd_slp, rho*sd_int*sd_slp, sd_slp^2), 2))
subj$b0 <- re[, 1]; subj$b1 <- re[, 2]

long <- subj %>% tidyr::crossing(AVISITN = weeks) %>%
  mutate(YEARS = AVISITN / 52.18, CENT10 = AMYLCENT / 10,
         mu = beta["(Intercept)"] + beta["CENT10"]*CENT10 + beta["YEARS"]*YEARS + beta["CENT10:YEARS"]*CENT10*YEARS +
              beta["AGEYR"]*AGEYR + beta["EDCCNTU"]*EDCCNTU + beta[grep("APOE4", names(beta))]*(APOE4 == "carrier") +
              # a small simulated treatment effect so the MMRM tab has something to show
              if_else(COHORT == "Active", 0.12 * YEARS, 0),
         AVAL = round(mu + b0 + b1*YEARS + rnorm(n(), 0, sd_res), 3))

# MAR-type dropout: hazard increases with age and with low last-observed PACC
drop <- long %>% arrange(BID, AVISITN) %>% group_by(BID) %>%
  mutate(p_drop = plogis(-4 + 0.03*(AGEYR - 72) - 0.15*lag(AVAL, default = first(AVAL))),
         u = runif(n()), dropped = cumsum(u < p_drop & AVISITN > 0) > 0) %>% ungroup() %>% filter(!dropped)

adqs <- drop %>% group_by(BID) %>% mutate(BASE = AVAL[AVISITN == 0][1], CHG = AVAL - BASE) %>% ungroup() %>%
  mutate(AVISIT = factor(paste0("Week ", AVISITN), levels = paste0("Week ", weeks)),
         COHORT = factor(COHORT, levels = names(n))) %>%
  select(BID, SITE, COHORT, AVISIT, AVISITN, AVAL, BASE, CHG, AGEYR, EDCCNTU, SEX, APOE4, AMYLCENT)
saveRDS(adqs, "shiny/data/synthetic_adqs.rds")
message("Synthetic ADQS: ", nrow(adqs), " rows, ", n_distinct(adqs$BID), " subjects (no real data)")

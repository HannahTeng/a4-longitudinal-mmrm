# ------------------------------------------------------------------
# 05_tables_figures.R: Table 1, observed and LS-mean figures, trajectories by amyloid group
# ------------------------------------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_minimal(base_size = 12))

adsl <- readRDS(file.path(DERIVED_DIR, "adsl.rds"))
adqs <- readRDS(file.path(DERIVED_DIR, "adqs_pacc.rds"))

# ---- Table 1: baseline characteristics, A4 mITT by arm + LEARN
t1 <- adsl %>%
  filter((SUBSTUDY == "A4" & MITTFL == 1) | SUBSTUDY == "LEARN") %>%
  mutate(grp = if_else(SUBSTUDY == "LEARN", "LEARN", as.character(TX))) %>%
  group_by(grp) %>%
  summarise(N = as.character(n()),
            `Age, mean (SD)`          = sprintf("%.1f (%.1f)", mean(AGEYR, na.rm=TRUE), sd(AGEYR, na.rm=TRUE)),
            `Female, n (%)`           = sprintf("%d (%.1f)", sum(SEX == "Female", na.rm=TRUE), 100*mean(SEX == "Female", na.rm=TRUE)),
            `Education yr, mean (SD)` = sprintf("%.1f (%.1f)", mean(EDCCNTU, na.rm=TRUE), sd(EDCCNTU, na.rm=TRUE)),
            `APOE4 carrier, n (%)`    = sprintf("%d (%.1f)", sum(APOE4 == "carrier", na.rm=TRUE), 100*mean(APOE4 == "carrier", na.rm=TRUE)),
            `Amyloid SUVR, mean (SD)` = sprintf("%.2f (%.2f)", mean(SUVRCER, na.rm=TRUE), sd(SUVRCER, na.rm=TRUE)),
            `Centiloids, mean (SD)`   = sprintf("%.1f (%.1f)", mean(AMYLCENT, na.rm=TRUE), sd(AMYLCENT, na.rm=TRUE)),
            `MMSE, mean (SD)`         = sprintf("%.1f (%.1f)", mean(MMSEBL, na.rm=TRUE), sd(MMSEBL, na.rm=TRUE))) %>%
  pivot_longer(-grp, names_to = "Characteristic", values_to = "v") %>%
  pivot_wider(names_from = grp, values_from = v) %>%
  select(Characteristic, Placebo, Solanezumab, LEARN)
write_csv(t1, file.path(OUT_TAB, "table1_baseline.csv"))

# ---- Figure 1: observed mean change from baseline (+/- SE) by arm, double-blind period
f1 <- adqs %>% filter(SUBSTUDY == "A4", MITTFL == 1, !is.na(CHG)) %>%
  group_by(TX, AVISITN) %>% summarise(m = mean(CHG), se = sd(CHG)/sqrt(n()), n = n(), .groups = "drop")
p1 <- ggplot(f1, aes(AVISITN, m, colour = TX)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 4) + geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = PACC_VISIT_WEEKS) +
  labs(x = "Week", y = "Mean change from baseline in PACC (± SE)", colour = NULL,
       title = "A4: observed PACC change from baseline, mITT, double-blind period") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "fig1_observed_change_by_arm.png"), p1, width = 8, height = 5, dpi = 150)

# ---- Figure 2: MMRM LS-means by arm
lsm <- read_csv(file.path(OUT_TAB, "aim1_mmrm_lsmeans.csv"), show_col_types = FALSE) %>%
  mutate(week = as.integer(gsub("Week ", "", AVISIT)))
p2 <- ggplot(lsm, aes(week, emmean, colour = TX)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 4) + geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = PACC_VISIT_WEEKS) +
  labs(x = "Week", y = "LS-mean change from baseline in PACC (95% CI)", colour = NULL,
       title = "A4 primary analysis (MMRM, unstructured covariance, Kenward-Roger df)") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "fig2_mmrm_lsmeans_by_arm.png"), p2, width = 8, height = 5, dpi = 150)

# ---- Figure 3: spaghetti + loess by amyloid group (placebo + LEARN)
obs <- adqs %>% filter(COHORT %in% c("A4 Placebo", "LEARN (Abeta-)"), !is.na(AVAL), !is.na(AMYLCENT)) %>%
  mutate(AMYL_GRP = cut(AMYLCENT, c(-Inf, 20, 50, Inf), right = FALSE, labels = c("<20 CL", "20 to <50 CL", ">=50 CL")))
# Group-level smooths only: individual profiles are row-level DUA data and are not plotted in committed figures.
p3 <- ggplot(obs, aes(AVISITN / 52.18, AVAL)) +
  geom_smooth(aes(colour = AMYL_GRP), method = "loess", se = TRUE, span = 1) +
  labs(x = "Years from baseline", y = "PACC composite score", colour = "Baseline amyloid",
       title = "PACC trajectories by baseline amyloid, untreated participants",
       subtitle = "A4 placebo and LEARN; loess smooth with 95% band") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "fig3_trajectories_by_amyloid.png"), p3, width = 8, height = 5.5, dpi = 150)

# ---- Figure 4: Aim 2 LS-means by amyloid group
lsg <- read_csv(file.path(OUT_TAB, "aim2_mmrm_lsmeans_by_amyloid.csv"), show_col_types = FALSE) %>%
  mutate(week = as.integer(gsub("Week ", "", AVISIT)))
p4 <- ggplot(lsg, aes(week, emmean, colour = AMYL_GRP)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 4) + geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = PACC_VISIT_WEEKS) +
  labs(x = "Week", y = "LS-mean change from baseline in PACC (95% CI)", colour = "Baseline amyloid",
       title = "MMRM: PACC change by baseline amyloid group (untreated participants)") +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "fig4_mmrm_lsmeans_by_amyloid.png"), p4, width = 8, height = 5, dpi = 150)
message("Tables and figures written.")

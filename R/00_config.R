# ------------------------------------------------------------------
# 00_config.R: paths, population flags, analysis constants
# A4 Longitudinal Cognitive Decline: MMRM analysis pipeline
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(forcats)
})

RAW_DIR     <- "data/raw"       # drop the Synapse CSVs here (never committed)
DERIVED_DIR <- "data/derived"
OUT_TAB     <- "output/tables"
OUT_FIG     <- "output/figures"

# Visit schedule used in the MMRM (PACC administered every 24 weeks in A4)
PACC_VISIT_WEEKS <- c(0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240)
PRIMARY_WEEK     <- 240   # A4 primary endpoint: change from baseline in PACC at week 240

# Analysis-population definitions (mirrors SUBJINFO flags)
#   mITT  : randomized A4 participants with >= 1 post-baseline PACC (MITTFL == 1)
#   OBS   : A4 placebo arm + LEARN cohort (for amyloid x time modeling)
POP_MITT <- function(adsl) adsl %>% filter(SUBSTUDY == "A4", MITTFL == 1)
POP_OBS  <- function(adsl) adsl %>% filter((SUBSTUDY == "A4" & TX == "Placebo") | SUBSTUDY == "LEARN")

read_raw <- function(name) read_csv(file.path(RAW_DIR, paste0(name, ".csv")),
                                    show_col_types = FALSE, guess_max = 1e5)

# A4 Longitudinal Cognitive Decline: MMRM Analysis Pipeline

Reproducible analysis of longitudinal cognition (Preclinical Alzheimer Cognitive Composite, PACC) in the A4 Study
(Anti-Amyloid Treatment in Asymptomatic Alzheimer's Disease) and its companion LEARN cohort. The primary analysis
is a mixed model for repeated measures (MMRM). It is implemented twice, in R and in SAS, as independent
programs (double programming), and the two sets of results are reconciled numerically. The project is also the
Final Data Analysis Project for UCLA BIOSTAT M236 (Longitudinal Data, Fall 2026).

> Data are not included. A4/LEARN data are distributed under a data use agreement via the ATRI data
> sharing platform on Synapse (`syn61250768`). Place the public-release CSVs in `data/raw/` and the pipeline
> runs end-to-end. Only code, logs and aggregate outputs (tables, figures) are committed.

## Two-engine design (R and SAS)

| Engine | Aim 1 program | Role |
|---|---|---|
| R | `R/03_mmrm_primary.R` (`mmrm::mmrm`; `nlme::gls` as an internal cross-check) | Primary implementation, derivations, figures, dashboard |
| SAS | `sas/01_mmrm_primary.sas` (`PROC MIXED`, `REPEATED / TYPE=UN`, `DDFM=KR`) | Independent re-implementation from the exported analysis datasets |

The SAS program does not read any R model output; it reads only the exported ADQS data and refits the model.
`R/07_reconcile_r_vs_sas.R` then compares the visit-specific treatment-difference estimate, standard error and
Kenward-Roger denominator df. For Aim 1 the maximum absolute differences across the ten visits are 1.2 x 10^-5
(estimate), 9 x 10^-4 (SE) and 0.14 (df); the df gap reflects small differences in the Kenward-Roger
implementations and does not change any estimate or conclusion. The SAS program for Aim 2
(`sas/02_amyloid_slope.sas`) is provided; its R-versus-SAS reconciliation has not yet been run.

## Aims

| Aim | Question | Population | Model |
|---|---|---|---|
| 1 | Does solanezumab slow PACC decline vs placebo? (replication of the published primary analysis) | A4 mITT, double-blind period, weeks 24–240 | MMRM on change from baseline: `CHG ~ TX*AVISIT + BASE + AGE + EDU + APOE4`, unstructured covariance, Kenward-Roger df |
| 2 | How does baseline amyloid burden shape the PACC trajectory in untreated participants? | A4 placebo arm + LEARN (amyloid-negative) | (a) random intercept/slope LMM `PACC ~ centiloids × years`; (b) MMRM by amyloid group; (c) covariance-structure selection (CS / AR(1) / ARH(1) / UN) by REML AIC/BIC and LRT |

## Results (A4/LEARN public release v1.2)

Aim 1, primary endpoint replication. Solanezumab minus placebo difference in change from baseline PACC at
week 240: -0.32 (95% CI -0.93 to 0.29), p = 0.30, n = 1,115. The published A4 result (Sperling et al.,
*NEJM* 2023) was -0.30 (-0.91 to 0.31), p = 0.33, n = 1,147. The small n difference is because 31 mITT
participants have post-baseline PACC only at an early-termination visit (code 999), which the public
release does not map to a scheduled week; the released analysis windowed those into scheduled visits.
Two independent implementations (`mmrm::mmrm` and `nlme::gls` with `corSymm` + `varIdent`) agree to
< 2 × 10⁻⁴ on every fixed-effect coefficient.

Aim 2, amyloid and decline. In untreated participants (n = 1,127), each additional 10 centiloids of
baseline amyloid steepens the annual PACC slope by -0.058 (SE 0.005, t = -11.1). The unstructured
covariance fits far better than CS, AR(1) or ARH(1) (ΔAIC > 1,500; LRT vs ARH(1) p < 10⁻³⁰⁰), which is the
empirical justification for the unstructured MMRM used in the primary analysis.

Figures: `output/figures/` (observed means by arm, MMRM LS-means by arm, individual trajectories by amyloid
group, MMRM LS-means by amyloid group). Tables: `output/tables/`.

## Pipeline

```
R/00_config.R              paths, visit schedule, population definitions
R/01_derive_adsl.R         SUBJINFO -> ADSL (subject-level; flags, baseline covariates)   + derivation QC
R/02_derive_adqs.R         PACC.csv -> ADQS (one row per subject x scheduled visit)      + QC vs released CHG
R/03_mmrm_primary.R        Aim 1 MMRM (mmrm) + independent gls cross-check
R/04_amyloid_slope.R       Aim 2 LMM, MMRM by amyloid group, covariance comparison
R/05_tables_figures.R      Table 1, figures
R/06_export_for_sas.R      CSV exports for SAS
R/07_reconcile_r_vs_sas.R  R vs SAS reconciliation report (reads the SAS CSV outputs)
sas/00_import.sas          import the exported CSVs into SAS datasets
sas/01_mmrm_primary.sas    PROC MIXED re-implementation of Aim 1 (TYPE=UN, DDFM=KR)
sas/02_amyloid_slope.sas   PROC MIXED re-implementation of Aim 2 (random slopes; CS/AR(1)/ARH(1)/UN)
run_all.R                  runs R/01–05
```

```r
# R >= 4.3; packages: dplyr tidyr readr stringr forcats ggplot2 nlme lme4 emmeans mmrm
source("run_all.R")
source("R/06_export_for_sas.R")   # then run sas/*.sas (e.g. on SAS OnDemand for Academics)
source("R/07_reconcile_r_vs_sas.R")
```

## Data model

* ADSL (`data/derived/adsl.rds`): one row per screened participant. `SUBSTUDY` (A4 / LEARN / SF),
  `TX`, `MITTFL`, age, education, sex, APOE4 (`APOE4_AN` = analysis flag, missing → non-carrier per the
  randomization plan), florbetapir SUVR and centiloids (`AMYLCENT = 183.07 × SUVR − 177.26`).
* ADQS (`data/derived/adqs_pacc.rds`): PACC at baseline (visit 6) and every 24 weeks through week 240
  (visit 66), `AVAL` = PACC composite, `BASE`, `CHG`; open-label (visits > 66) and unscheduled / early-termination
  visits (701, 997, 999) are excluded from the double-blind analyses.
* QC gates (`stopifnot`): every randomized participant has an arm; mITT ⊂ randomized; recomputed `CHG`
  matches the released change-from-baseline to 10⁻⁵; no duplicate subject-visit rows.

## Methods notes

* PACC = mean of z-scores of FCSRT free+total recall, Logical Memory delayed recall, Digit Symbol, and
  MMSE, standardized to the screening population (Donohue et al., *JAMA Neurol* 2014).
* MMRM treats visit as categorical, estimates a separate mean at each visit, and uses REML with an
  unstructured covariance so no parametric time trend or correlation pattern is assumed; missing data are
  handled under MAR. Kenward-Roger degrees of freedom match the SAS `DDFM=KR` default used in trial reporting.
* The random-slope model in Aim 2 is the continuous-time analogue (centiloids as a continuous moderator of
  slope); the MMRM-by-tertile analysis is its categorical-time counterpart.

## References

* Sperling RA, et al. Trial of Solanezumab in Preclinical Alzheimer's Disease. *N Engl J Med* 2023;389:1096-1107.
* Donohue MC, et al. The Preclinical Alzheimer Cognitive Composite. *JAMA Neurol* 2014;71:961-970.
* A4/LEARN data: https://www.synapse.org/Synapse:syn61250768 (data use agreement required).

## Shiny dashboard

`shiny/app.R` is a longitudinal monitoring and analysis dashboard. Public deployment:
https://hannahteng.shinyapps.io/pacc-longitudinal-monitor/

| Tab | Content |
|---|---|
| Overview | Purpose, model specification, description of the R/SAS double-programming, upload instructions |
| Enrollment & visits | Participants, retention, visit completion by cohort, baseline characteristics |
| Trajectories | Observed PACC profiles and group means; mean change from baseline with SE; stratification by cohort, APOE4, sex or baseline amyloid |
| MMRM | Live fit (`mmrm`, unstructured covariance, Kenward-Roger df); LS-means and cohort contrasts by visit with 95% CI |
| Covariance | CS, AR(1), ARH(1) and unstructured compared by REML AIC/BIC; LRT of ARH(1) vs unstructured |
| Validation | R (`mmrm`) vs SAS (`PROC MIXED`) estimate, SE and df by visit, with absolute differences |
| Data | Filtered analysis dataset |

Data sources, in order of preference: an uploaded CSV; the local derived A4/LEARN data (`data/derived/adqs_pacc.rds`,
when present); `shiny/data/synthetic_adqs.rds`. The public deployment serves only the synthetic dataset, simulated
from the fitted Aim 2 model (`shiny/make_synthetic.R`) with the structure of the study and no real participant data.
Run locally with `shiny::runApp("shiny")`; set `PACC_DASH_SYNTHETIC=1` to force the synthetic data.

Upload format: long CSV, one row per subject and visit. Required columns `BID, COHORT, AVISITN, AVAL, AGEYR, EDCCNTU,
SEX, APOE4`; optional `SITE, BASE, CHG, AMYLCENT` (`BASE` and `CHG` are derived from the week 0 value when absent).
`SEX` is `Female`/`Male`, `APOE4` is `carrier`/`non-carrier`, a baseline visit coded 0 is required, and the file is
limited to 5 MB, 60,000 rows, 15 visits and 6 cohorts. Uploaded files are held in memory for the session only; do not
upload identifiable data or data covered by a data use agreement to the hosted app.

Deployment: `source("shiny/deploy.R")` after registering shinyapps.io credentials with `rsconnect::setAccountInfo()`.
The deployment bundle contains `app.R`, the synthetic dataset and the aggregate reconciliation table only.

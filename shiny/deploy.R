# One-time: rsconnect::setAccountInfo(name = "hannahteng", token = "<token>", secret = "<secret>")
#           (Account -> Tokens on shinyapps.io). Then from the repo root:
#   source("shiny/deploy.R")
# Only app.R, the synthetic dataset and the aggregate R-vs-SAS reconciliation table are bundled; ../data is NOT uploaded,
# so the public app always runs on simulated data.
rsconnect::deployApp(appDir = "shiny", appName = "pacc-longitudinal-monitor",
                     appFiles = c("app.R", "data/synthetic_adqs.rds", "data/aim1_reconciliation_r_vs_sas.csv"), forceUpdate = TRUE)

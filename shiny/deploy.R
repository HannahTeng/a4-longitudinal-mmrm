# One-time: rsconnect::setAccountInfo(name = "hannahteng", token = "<token>", secret = "<secret>")
#           (Account -> Tokens on shinyapps.io). Then from the repo root:
#   source("shiny/deploy.R")
# Only shiny/app.R and shiny/data/synthetic_adqs.rds are bundled; ../data is NOT uploaded,
# so the public app always runs on simulated data.
rsconnect::deployApp(appDir = "shiny", appName = "pacc-longitudinal-monitor",
                     appFiles = c("app.R", "data/synthetic_adqs.rds"), forceUpdate = TRUE)

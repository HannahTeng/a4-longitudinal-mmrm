# ------------------------------------------------------------------
# Longitudinal Trial Monitoring Dashboard — PACC in a preclinical AD study
# Hannah (Zihan) Teng · github.com/hannahteng/a4-longitudinal-mmrm
#
# Data source (in order of preference):
#   1. ../data/derived/adqs_pacc.rds  — real A4/LEARN derived data (local use only, never deployed)
#   2. data/synthetic_adqs.rds        — simulated data with A4-like structure (what the public app serves)
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(dplyr); library(tidyr); library(ggplot2)
  library(plotly); library(DT); library(mmrm); library(emmeans); library(nlme)
})

# ---------- data ----------
load_data <- function() {
  real <- file.path("..", "data", "derived", "adqs_pacc.rds")
  if (file.exists(real) && Sys.getenv("PACC_DASH_SYNTHETIC") != "1") {
    d <- readRDS(real) %>%
      filter(COHORT != "Screen failure") %>%
      transmute(BID, SITE = NA_character_, COHORT = droplevels(factor(COHORT)), AVISIT, AVISITN,
                AVAL, BASE, CHG, AGEYR, EDCCNTU, SEX, APOE4 = APOE4_AN, AMYLCENT)
    list(data = d, label = "A4/LEARN derived data (local, not for distribution)", real = TRUE)
  } else {
    list(data = readRDS(file.path("data", "synthetic_adqs.rds")),
         label = "Simulated data with A4-like structure (no real participants)", real = FALSE)
  }
}
SRC   <- load_data()
ADQS  <- SRC$data
WEEKS <- sort(unique(ADQS$AVISITN))
pal   <- c("#2F6DB5", "#E07A2D", "#4F9A6B", "#8C5AA8")

# ---------- UI ----------
ui <- page_navbar(
  title = "PACC Longitudinal Monitoring",
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_google("Inter")),
  sidebar = sidebar(
    width = 290,
    tags$small(class = "text-muted", SRC$label),
    checkboxGroupInput("cohorts", "Cohorts", choices = levels(ADQS$COHORT), selected = levels(ADQS$COHORT)),
    sliderInput("maxweek", "Follow-up window (weeks)", min = 24, max = max(WEEKS), value = max(WEEKS), step = 24),
    selectInput("strat", "Stratify trajectories by",
                c("Cohort" = "COHORT", "APOE4" = "APOE4", "Baseline amyloid group" = "AMYL_GRP", "Sex" = "SEX")),
    hr(),
    tags$small("MMRM: change from baseline ~ cohort × visit + baseline + age + education + APOE4, ",
               "unstructured covariance, Kenward-Roger df (R package ", tags$code("mmrm"), ").")
  ),

  nav_panel("Enrollment & visits",
    layout_columns(
      fill = FALSE,
      value_box("Participants", textOutput("n_subj"), theme = "primary"),
      value_box("PACC assessments", textOutput("n_obs"), theme = "secondary"),
      value_box(paste0("Retention at week ", max(WEEKS)), textOutput("retention"), theme = "success"),
      value_box("Median follow-up (wk)", textOutput("med_fu"), theme = "info")
    ),
    layout_columns(
      card(card_header("Visit completion by cohort (% of baseline)"), plotlyOutput("p_completion", height = 340)),
      card(card_header("Baseline characteristics"), DTOutput("t_baseline"))
    )
  ),

  nav_panel("Trajectories",
    layout_columns(
      card(card_header("Individual profiles and group means (observed PACC)"), plotlyOutput("p_spaghetti", height = 460)),
      card(card_header("Observed mean change from baseline (± SE)"), plotlyOutput("p_meanchg", height = 460))
    )
  ),

  nav_panel("MMRM",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Fit"),
           selectInput("ref", "Reference cohort", levels(ADQS$COHORT)),
           actionButton("fit", "Fit MMRM", class = "btn-primary"),
           br(), br(), verbatimTextOutput("fit_msg")),
      card(card_header("LS-mean change from baseline by cohort (95% CI)"), plotlyOutput("p_lsm", height = 400))
    ),
    card(card_header("Cohort differences vs reference, by visit"), DTOutput("t_diff"))
  ),

  nav_panel("Covariance",
    card(card_header("Covariance-structure comparison for the marginal model (REML; lower AIC is better)"),
         p(class = "text-muted", "Fits the same mean model with compound symmetry, AR(1), heterogeneous AR(1) and unstructured ",
           "covariance via nlme::gls and reports AIC/BIC and the LRT of ARH(1) vs unstructured."),
         actionButton("fitcov", "Run comparison", class = "btn-primary"), br(), br(),
         DTOutput("t_cov"), verbatimTextOutput("cov_lrt"))
  ),

  nav_panel("Data", card(DTOutput("t_data")))
)

# ---------- server ----------
server <- function(input, output, session) {

  dat <- reactive({
    ADQS %>% filter(COHORT %in% input$cohorts, AVISITN <= input$maxweek) %>%
      mutate(AMYL_GRP = cut(AMYLCENT, c(-Inf, 20, 50, Inf), labels = c("<20 CL", "20-50 CL", ">50 CL")),
             AVISIT = droplevels(AVISIT), COHORT = droplevels(COHORT))
  })

  # --- KPIs
  output$n_subj <- renderText(format(n_distinct(dat()$BID), big.mark = ","))
  output$n_obs  <- renderText(format(nrow(dat()), big.mark = ","))
  output$retention <- renderText({
    d <- dat(); last <- max(d$AVISITN)
    sprintf("%.0f%%", 100 * n_distinct(d$BID[d$AVISITN == last]) / n_distinct(d$BID))
  })
  output$med_fu <- renderText({ dat() %>% group_by(BID) %>% summarise(m = max(AVISITN)) %>% pull(m) %>% median() })

  output$p_completion <- renderPlotly({
    d <- dat() %>% count(COHORT, AVISITN) %>% group_by(COHORT) %>% mutate(pct = 100 * n / n[AVISITN == 0]) %>% ungroup()
    g <- ggplot(d, aes(AVISITN, pct, colour = COHORT, text = paste0(COHORT, "<br>Week ", AVISITN, ": n = ", n))) +
      geom_line(aes(group = COHORT)) + geom_point() + scale_colour_manual(values = pal) +
      scale_x_continuous(breaks = WEEKS) + labs(x = "Week", y = "% of baseline participants", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  output$t_baseline <- renderDT({
    dat() %>% filter(AVISITN == 0) %>% group_by(Cohort = COHORT) %>%
      summarise(N = n(), `Age` = sprintf("%.1f (%.1f)", mean(AGEYR), sd(AGEYR)),
                `Education` = sprintf("%.1f (%.1f)", mean(EDCCNTU, na.rm = TRUE), sd(EDCCNTU, na.rm = TRUE)),
                `Female %` = sprintf("%.0f", 100 * mean(SEX == "Female", na.rm = TRUE)),
                `APOE4 carrier %` = sprintf("%.0f", 100 * mean(APOE4 == "carrier", na.rm = TRUE)),
                `Centiloids` = sprintf("%.1f (%.1f)", mean(AMYLCENT, na.rm = TRUE), sd(AMYLCENT, na.rm = TRUE)),
                `Baseline PACC` = sprintf("%.2f (%.2f)", mean(BASE, na.rm = TRUE), sd(BASE, na.rm = TRUE))) %>%
      datatable(options = list(dom = "t"), rownames = FALSE, caption = "Mean (SD) unless noted")
  })

  # --- trajectories
  output$p_spaghetti <- renderPlotly({
    d <- dat() %>% filter(!is.na(.data[[input$strat]]))
    set.seed(7); samp <- sample(unique(d$BID), min(120, n_distinct(d$BID)))
    m <- d %>% group_by(grp = .data[[input$strat]], AVISITN) %>% summarise(m = mean(AVAL), .groups = "drop")
    g <- ggplot() +
      geom_line(data = filter(d, BID %in% samp), aes(AVISITN, AVAL, group = BID), colour = "grey75", alpha = 0.5) +
      geom_line(data = m, aes(AVISITN, m, colour = grp), linewidth = 1.2) +
      geom_point(data = m, aes(AVISITN, m, colour = grp, text = sprintf("%s<br>Week %d: mean PACC %.2f", grp, AVISITN, m))) +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = WEEKS) +
      labs(x = "Week", y = "PACC composite", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  output$p_meanchg <- renderPlotly({
    d <- dat() %>% filter(AVISITN > 0, !is.na(CHG), !is.na(.data[[input$strat]])) %>%
      group_by(grp = .data[[input$strat]], AVISITN) %>% summarise(m = mean(CHG), se = sd(CHG) / sqrt(n()), n = n(), .groups = "drop")
    g <- ggplot(d, aes(AVISITN, m, colour = grp, text = sprintf("%s<br>Week %d: %.2f ± %.2f (n=%d)", grp, AVISITN, m, se, n))) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 5) + geom_line(aes(group = grp)) + geom_point() +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = WEEKS) +
      labs(x = "Week", y = "Mean change from baseline", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  # --- MMRM
  fit <- eventReactive(input$fit, {
    d <- dat() %>% filter(AVISITN > 0, !is.na(CHG), !is.na(BASE), !is.na(AGEYR), !is.na(EDCCNTU), !is.na(APOE4)) %>%
      mutate(BID = factor(BID), COHORT = relevel(COHORT, input$ref))
    validate(need(nlevels(d$COHORT) >= 2, "Select at least two cohorts."))
    withProgress(message = "Fitting MMRM (unstructured covariance)...", {
      f <- mmrm(CHG ~ COHORT * AVISIT + BASE + AGEYR + EDCCNTU + APOE4 + us(AVISIT | BID), data = d, method = "Kenward-Roger")
    })
    em <- emmeans(f, ~ COHORT | AVISIT)
    list(fit = f, n = n_distinct(d$BID), rows = nrow(d),
         lsm = as.data.frame(em) %>% mutate(week = as.integer(gsub("Week ", "", AVISIT))),
         diff = as.data.frame(contrast(em, "trt.vs.ctrl")) %>%
           mutate(lower = estimate - qt(0.975, df) * SE, upper = estimate + qt(0.975, df) * SE))
  })
  output$fit_msg <- renderPrint({
    f <- fit()
    cat(sprintf("Subjects: %d\nPost-baseline records: %d\nAIC: %.1f\n",
                f$n, f$rows, AIC(f$fit)))
  })
  output$p_lsm <- renderPlotly({
    l <- fit()$lsm
    g <- ggplot(l, aes(week, emmean, colour = COHORT, text = sprintf("%s<br>Week %d: %.2f (%.2f, %.2f)", COHORT, week, emmean, lower.CL, upper.CL))) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 5) + geom_line(aes(group = COHORT)) + geom_point(size = 2) +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = WEEKS) +
      labs(x = "Week", y = "LS-mean change from baseline", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })
  output$t_diff <- renderDT({
    fit()$diff %>% transmute(Contrast = contrast, Visit = AVISIT, Estimate = round(estimate, 3), SE = round(SE, 3),
                             df = round(df, 1), `95% CI` = sprintf("(%.3f, %.3f)", lower, upper), p = signif(p.value, 2)) %>%
      datatable(options = list(pageLength = 12, dom = "tp"), rownames = FALSE)
  })

  # --- covariance comparison
  covcmp <- eventReactive(input$fitcov, {
    d <- dat() %>% filter(AVISITN > 0, !is.na(CHG), !is.na(BASE), !is.na(AGEYR), !is.na(EDCCNTU), !is.na(APOE4)) %>%
      mutate(BID = factor(BID))
    mm <- CHG ~ COHORT * AVISIT + BASE + AGEYR + EDCCNTU + APOE4
    tv <- ~ as.integer(AVISIT) | BID
    withProgress(message = "Fitting four covariance structures...", {
      fits <- list(
        "Compound symmetry"    = gls(mm, d, correlation = corCompSymm(form = tv), method = "REML"),
        "AR(1)"                = gls(mm, d, correlation = corAR1(form = tv), method = "REML"),
        "Heterogeneous AR(1)"  = gls(mm, d, correlation = corAR1(form = tv), weights = varIdent(form = ~ 1 | AVISIT), method = "REML"),
        "Unstructured"         = gls(mm, d, correlation = corSymm(form = tv), weights = varIdent(form = ~ 1 | AVISIT),
                                     method = "REML", control = glsControl(maxIter = 200, msMaxIter = 200)))
    })
    tab <- tibble(Structure = names(fits), `Cov. parameters` = sapply(fits, function(f) length(unlist(f$modelStruct))),
                  AIC = round(sapply(fits, AIC), 1), BIC = round(sapply(fits, BIC), 1), logLik = round(sapply(fits, logLik), 1)) %>% arrange(AIC)
    list(tab = tab, lrt = anova(fits[["Heterogeneous AR(1)"]], fits[["Unstructured"]]))
  })
  output$t_cov <- renderDT(datatable(covcmp()$tab, options = list(dom = "t"), rownames = FALSE))
  output$cov_lrt <- renderPrint({ a <- covcmp()$lrt; cat(sprintf("LRT heterogeneous AR(1) vs unstructured: chi-sq = %.1f on %d df, p = %.3g\n",
                                                                  a$L.Ratio[2], a$df[2] - a$df[1], a$`p-value`[2])) })

  output$t_data <- renderDT(datatable(dat() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
                                      options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE, filter = "top"))
}

shinyApp(ui, server)

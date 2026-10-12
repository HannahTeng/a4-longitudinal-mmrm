# ------------------------------------------------------------------
# PACC Longitudinal Monitoring: MMRM dashboard for a preclinical AD study
# Source: github.com/hannahteng/a4-longitudinal-mmrm
#
# Data source, in order of preference:
#   1. A user-uploaded CSV (session-local, never persisted)
#   2. ../data/derived/adqs_pacc.rds, the real A4/LEARN derived data (local use only, never deployed)
#   3. data/synthetic_adqs.rds, simulated data with A4-like structure (served by the public app)
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(dplyr); library(tidyr); library(ggplot2)
  library(plotly); library(DT); library(mmrm); library(emmeans); library(nlme)
})

options(shiny.maxRequestSize = 5 * 1024^2)

# ---------- data ----------
load_default <- function() {
  real <- file.path("..", "data", "derived", "adqs_pacc.rds")
  if (file.exists(real) && Sys.getenv("PACC_DASH_SYNTHETIC") != "1") {
    d <- readRDS(real) %>%
      filter(COHORT != "Screen failure") %>%
      transmute(BID, SITE = NA_character_, COHORT = droplevels(factor(COHORT)), AVISIT, AVISITN,
                AVAL, BASE, CHG, AGEYR, EDCCNTU, SEX, APOE4 = APOE4_AN, AMYLCENT)
    list(data = d, label = "A4/LEARN derived data (local use, not for distribution)", real = TRUE)
  } else {
    list(data = readRDS(file.path("data", "synthetic_adqs.rds")),
         label = "Simulated data with A4-like structure (no real participants)", real = FALSE)
  }
}
DEFAULT   <- load_default()
SYNTHETIC <- readRDS(file.path("data", "synthetic_adqs.rds"))   # template rows always come from simulated data

REQUIRED_COLS <- c("BID", "COHORT", "AVISITN", "AVAL", "AGEYR", "EDCCNTU", "SEX", "APOE4")
OPTIONAL_COLS <- c("SITE", "BASE", "CHG", "AMYLCENT")

# Validate and standardize an uploaded long-format dataset (one row per subject x visit).
# Returns list(data = <tibble>) on success or list(error = <character>) on failure.
prepare_upload <- function(raw) {
  fail <- function(...) list(error = paste0(...))
  names(raw) <- toupper(trimws(names(raw)))
  miss <- setdiff(REQUIRED_COLS, names(raw))
  if (length(miss)) return(fail("Missing required column(s): ", paste(miss, collapse = ", "), "."))
  if (nrow(raw) > 60000) return(fail("File has ", nrow(raw), " rows; the limit is 60,000."))

  d <- raw[, intersect(c(REQUIRED_COLS, OPTIONAL_COLS), names(raw)), drop = FALSE]
  for (v in intersect(c("AVISITN", "AVAL", "AGEYR", "EDCCNTU", "BASE", "CHG", "AMYLCENT"), names(d))) {
    if (all(is.na(d[[v]])) && v %in% REQUIRED_COLS) return(fail("Column ", v, " is entirely missing."))
    if (!is.numeric(d[[v]]) && !all(is.na(d[[v]]))) return(fail("Column ", v, " must be numeric."))
  }
  d$BID <- as.character(d$BID); d$COHORT <- as.character(d$COHORT)
  d$SEX <- as.character(d$SEX); d$APOE4 <- as.character(d$APOE4)

  if (anyNA(d$BID) || anyNA(d$COHORT) || anyNA(d$AVISITN)) return(fail("BID, COHORT and AVISITN may not be missing."))
  if (!all(na.omit(d$SEX) %in% c("Female", "Male"))) return(fail("SEX must be 'Female' or 'Male'."))
  if (!all(na.omit(d$APOE4) %in% c("carrier", "non-carrier"))) return(fail("APOE4 must be 'carrier' or 'non-carrier'."))
  if (anyDuplicated(d[c("BID", "AVISITN")])) return(fail("Duplicate subject-visit rows (BID x AVISITN) found."))
  if (any(tapply(d$COHORT, d$BID, dplyr::n_distinct) > 1)) return(fail("COHORT must be constant within subject."))
  nco <- dplyr::n_distinct(d$COHORT)
  if (nco < 2 || nco > 6) return(fail("COHORT must have between 2 and 6 levels (found ", nco, ")."))
  if (any(d$AVISITN < 0) || any(d$AVISITN != round(d$AVISITN))) return(fail("AVISITN must be a non-negative whole number of weeks."))
  wks <- sort(unique(d$AVISITN))
  if (!0 %in% wks) return(fail("AVISITN must include a baseline visit coded 0."))
  if (length(wks) < 2) return(fail("At least one post-baseline visit is required."))
  no_bl <- setdiff(unique(d$COHORT), unique(d$COHORT[d$AVISITN == 0]))
  if (length(no_bl)) return(fail("No baseline (week 0) records for cohort(s): ", paste(no_bl, collapse = ", "), "."))
  for (v in c("SEX", "APOE4")) if (dplyr::n_distinct(na.omit(d[[v]])) < 2) return(fail(v, " must take both values; a constant covariate cannot be adjusted for."))
  if (dplyr::n_distinct(na.omit(d$EDCCNTU)) < 2 || dplyr::n_distinct(na.omit(d$AGEYR)) < 2) return(fail("AGEYR and EDCCNTU must vary across subjects."))
  if (length(wks) > 15) return(fail("At most 15 distinct visits are supported (found ", length(wks), ")."))

  if (!"SITE" %in% names(d)) d$SITE <- NA_character_
  if (!"AMYLCENT" %in% names(d)) d$AMYLCENT <- NA_real_
  if (!"BASE" %in% names(d)) d <- d %>% group_by(BID) %>% mutate(BASE = AVAL[AVISITN == 0][1]) %>% ungroup()
  if (!"CHG" %in% names(d)) d$CHG <- d$AVAL - d$BASE
  d <- d %>%
    mutate(COHORT = factor(COHORT, levels = unique(COHORT)),
           AVISIT = factor(paste("Week", AVISITN), levels = paste("Week", wks))) %>%
    select(BID, SITE, COHORT, AVISIT, AVISITN, AVAL, BASE, CHG, AGEYR, EDCCNTU, SEX, APOE4, AMYLCENT)
  list(data = d)
}

# Reconciliation of R (mmrm) and SAS (PROC MIXED) treatment-difference estimates, bundled with the app.
RECON <- tryCatch(read.csv(file.path("data", "aim1_reconciliation_r_vs_sas.csv"), stringsAsFactors = FALSE),
                  error = function(e) NULL)

pal <- c("#2F6DB5", "#E07A2D", "#4F9A6B", "#8C5AA8", "#B5483B", "#6B7B8C")

# ---------- UI ----------
ui <- page_navbar(
  title = "PACC Longitudinal Monitoring",
  theme = bs_theme(version = 5, bootswatch = "flatly", base_font = font_google("Inter")),
  sidebar = sidebar(
    width = 300,
    textOutput("src_label") %>% tagAppendAttributes(class = "small text-muted"),
    checkboxGroupInput("cohorts", "Cohorts", choices = levels(DEFAULT$data$COHORT),
                       selected = levels(DEFAULT$data$COHORT)),
    selectInput("maxweek", "Follow-up through week",
                choices = sort(unique(DEFAULT$data$AVISITN))[-1], selected = max(DEFAULT$data$AVISITN)),
    selectInput("strat", "Stratify trajectories by",
                c("Cohort" = "COHORT", "APOE4" = "APOE4", "Baseline amyloid (centiloids)" = "AMYL_GRP", "Sex" = "SEX")),
    hr(),
    fileInput("upload", "Upload dataset (CSV, max 5 MB)", accept = ".csv", buttonLabel = "Browse", placeholder = "none"),
    uiOutput("upload_msg"),
    layout_columns(
      col_widths = c(6, 6),
      downloadButton("template", "Template", class = "btn-sm btn-outline-secondary w-100"),
      actionButton("reset", "Reset data", class = "btn-sm btn-outline-secondary w-100")
    ),
    hr(),
    tags$small("MMRM: change from baseline ~ cohort x visit + baseline + age + education + APOE4, ",
               "unstructured covariance, Kenward-Roger degrees of freedom (R package ", tags$code("mmrm"), ").")
  ),

  nav_panel("Overview",
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Purpose"),
        p("This dashboard supports monitoring and analysis of a longitudinal cognitive endpoint, the Preclinical Alzheimer ",
          "Cognitive Composite (PACC), in a randomized, parallel-group preclinical Alzheimer's disease study. It is built ",
          "around the analysis of the A4 Study and its companion LEARN cohort: the primary analysis is a mixed model for ",
          "repeated measures (MMRM) of change from baseline with visit as a categorical factor."),
        p("The public deployment runs on simulated data that mimic the structure of A4/LEARN (three cohorts, assessments every ",
          "24 weeks to week 240, dropout under a missing-at-random mechanism). It contains no real participants. A user-supplied ",
          "dataset can be analyzed instead (see below)."),
        h6("Statistical model"),
        tags$ul(
          tags$li("Response: change from baseline in PACC at each scheduled post-baseline visit."),
          tags$li("Fixed effects: cohort, visit (categorical), cohort x visit, baseline PACC, age, years of education, APOE4 carrier status."),
          tags$li("Covariance: unstructured within subject, estimated by REML; no parametric time trend or correlation pattern is assumed."),
          tags$li("Inference: Kenward-Roger degrees of freedom; LS-means and cohort contrasts by visit via ", tags$code("emmeans"), "."),
          tags$li("Missing data: ignorable (MAR) given the covariates and observed outcome history.")
        ),
        h6("Second engine: SAS"),
        p("The primary analysis model is fitted by two engines. R (", tags$code("mmrm::mmrm"),
          ") is the primary implementation; SAS (", tags$code("PROC MIXED"), " with ", tags$code("REPEATED / TYPE=UN"), " and ",
          tags$code("DDFM=KR"), ") independently refits the same model from the exported analysis dataset, without using any R model output. ",
          "Within R, ", tags$code("nlme::gls"), " with ", tags$code("corSymm"), " and ", tags$code("varIdent"),
          " is an additional check on the fixed effects. The R-versus-SAS comparison is on the Validation tab; the SAS programs are in ",
          "the repository under ", tags$code("sas/"), "."),
        p(class = "small text-muted", "Contrasts in the MMRM tab are unadjusted for multiplicity. LS-means weight categorical covariates ",
          "(APOE4) equally, as SAS LSMEANS does by default.")
      ),
      card(
        card_header("Tabs"),
        tags$dl(
          tags$dt("Enrollment & visits"), tags$dd("Participant counts, retention, visit completion by cohort, baseline characteristics."),
          tags$dt("Trajectories"), tags$dd("Observed PACC profiles and group means; mean change from baseline with standard errors."),
          tags$dt("MMRM"), tags$dd("Fit the model on the filtered data; LS-mean change by cohort and cohort contrasts by visit with 95% CI."),
          tags$dt("Covariance"), tags$dd("Compound symmetry, AR(1), heterogeneous AR(1) and unstructured structures compared by REML AIC/BIC, with an LRT of ARH(1) vs unstructured."),
          tags$dt("Validation"), tags$dd("R vs SAS reconciliation of the primary treatment-difference estimates."),
          tags$dt("Data"), tags$dd("Analysis dataset after filtering.")
        ),
        h6("Uploading data"),
        p("Upload a long-format CSV with one row per subject and visit. Required columns: ",
          tags$code(paste(REQUIRED_COLS, collapse = ", ")), ". Optional: ", tags$code(paste(OPTIONAL_COLS, collapse = ", ")),
          ". ", tags$code("BASE"), " and ", tags$code("CHG"), " are derived from the week 0 value when absent. ",
          "SEX is Female/Male; APOE4 is carrier/non-carrier. A baseline visit coded 0 and at most 15 visits and 6 cohorts are required. ",
          "Use the Template button for an example."),
        div(class = "alert alert-warning small",
            "Files are processed in memory for the current session and not stored, but the app runs on a third-party host. ",
            "Do not upload identifiable data or data restricted by a data use agreement (including A4/LEARN).")
      )
    )
  ),

  nav_panel("Enrollment & visits",
    layout_columns(
      fill = FALSE,
      value_box("Participants", textOutput("n_subj"), theme = "primary"),
      value_box("PACC assessments", textOutput("n_obs"), theme = "secondary"),
      value_box(textOutput("ret_title"), textOutput("retention"), theme = "success"),
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
      card(card_header("Observed mean change from baseline (+/- SE)"), plotlyOutput("p_meanchg", height = 460))
    )
  ),

  nav_panel("MMRM",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Fit"),
           selectInput("ref", "Reference cohort", levels(DEFAULT$data$COHORT)),
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

  nav_panel("Validation",
    card(card_header("Primary analysis: R (mmrm) vs SAS (PROC MIXED), Solanezumab - Placebo, change from baseline PACC"),
         p(class = "text-muted",
           "Each visit shows the treatment-difference estimate, standard error and Kenward-Roger denominator df from the two ",
           "independent implementations, and the absolute differences. The comparison was run on the A4 mITT analysis set; ",
           "only these aggregate estimates are shown."),
         DTOutput("t_recon"), br(), verbatimTextOutput("recon_summary"))
  ),

  nav_panel("Data", card(DTOutput("t_data")))
)

# ---------- server ----------
server <- function(input, output, session) {

  rv <- reactiveValues(src = DEFAULT, msg = NULL)
  adqs  <- reactive(rv$src$data)
  weeks <- reactive(sort(unique(adqs()$AVISITN)))

  fit_res <- reactiveVal(NULL)
  cov_res <- reactiveVal(NULL)

  # push a new dataset into the sidebar controls
  sync_inputs <- function(d) {
    lv <- levels(d$COHORT); wk <- sort(unique(d$AVISITN))
    updateCheckboxGroupInput(session, "cohorts", choices = lv, selected = lv)
    updateSelectInput(session, "ref", choices = lv, selected = lv[1])
    updateSelectInput(session, "maxweek", choices = wk[-1], selected = max(wk))
    ch <- c("Cohort" = "COHORT", "APOE4" = "APOE4", "Sex" = "SEX")
    if (any(!is.na(d$AMYLCENT))) ch <- c(ch, "Baseline amyloid (centiloids)" = "AMYL_GRP")
    updateSelectInput(session, "strat", choices = ch, selected = "COHORT")
    fit_res(NULL); cov_res(NULL)
  }

  observeEvent(input$upload, {
    res <- tryCatch(prepare_upload(read.csv(input$upload$datapath, stringsAsFactors = FALSE, check.names = FALSE)),
                    error = function(e) list(error = paste("Could not read file:", conditionMessage(e))))
    if (!is.null(res$error)) {
      rv$msg <- list(ok = FALSE, text = res$error)
    } else {
      rv$src <- list(data = res$data, label = paste0("Uploaded: ", input$upload$name), real = FALSE)
      rv$msg <- list(ok = TRUE, text = sprintf("Loaded %s rows, %s subjects, %d cohorts.",
                                               format(nrow(res$data), big.mark = ","),
                                               format(n_distinct(res$data$BID), big.mark = ","),
                                               n_distinct(res$data$COHORT)))
      sync_inputs(res$data)
    }
  })
  observeEvent(input$reset, {
    rv$src <- DEFAULT; rv$msg <- NULL; sync_inputs(DEFAULT$data)
  })
  observeEvent(input$cohorts, {
    sel <- if (isTRUE(input$ref %in% input$cohorts)) input$ref else input$cohorts[1]
    updateSelectInput(session, "ref", choices = input$cohorts, selected = sel)
  })
  output$upload_msg <- renderUI(if (!is.null(rv$msg)) div(class = paste("small", if (rv$msg$ok) "text-success" else "text-danger"), rv$msg$text))
  output$src_label  <- renderText(rv$src$label)

  output$template <- downloadHandler(
    filename = function() "pacc_upload_template.csv",
    content = function(file) {
      keep <- SYNTHETIC %>% distinct(BID, COHORT) %>% group_by(COHORT) %>% slice_head(n = 2) %>% pull(BID)
      d <- SYNTHETIC %>% filter(BID %in% keep) %>%
        select(BID, COHORT, AVISITN, AVAL, AGEYR, EDCCNTU, SEX, APOE4, SITE, AMYLCENT)
      write.csv(d, file, row.names = FALSE, na = "")
    })

  dat <- reactive({
    req(input$cohorts)
    adqs() %>% filter(COHORT %in% input$cohorts, AVISITN <= as.numeric(req(input$maxweek))) %>%
      mutate(AMYL_GRP = cut(AMYLCENT, c(-Inf, 20, 50, Inf), right = FALSE, labels = c("<20 CL", "20 to <50 CL", ">=50 CL")),
             AVISIT = droplevels(AVISIT), COHORT = droplevels(COHORT))
  })
  complete_cases <- reactive(
    dat() %>% filter(AVISITN > 0, !is.na(CHG), !is.na(BASE), !is.na(AGEYR), !is.na(EDCCNTU), !is.na(APOE4)) %>%
      mutate(BID = factor(BID), AVISIT = droplevels(AVISIT), COHORT = droplevels(COHORT)))

  # --- KPIs
  output$n_subj <- renderText(format(n_distinct(dat()$BID), big.mark = ","))
  output$n_obs  <- renderText(format(nrow(dat()), big.mark = ","))
  output$ret_title <- renderText(paste0("Retention at week ", max(dat()$AVISITN)))
  output$retention <- renderText({
    d <- dat(); last <- max(d$AVISITN)
    sprintf("%.0f%%", 100 * n_distinct(d$BID[d$AVISITN == last]) / n_distinct(d$BID))
  })
  output$med_fu <- renderText({ dat() %>% group_by(BID) %>% summarise(m = max(AVISITN)) %>% pull(m) %>% median() })

  output$p_completion <- renderPlotly({
    d <- dat() %>% count(COHORT, AVISITN) %>% group_by(COHORT) %>% mutate(pct = 100 * n / n[AVISITN == 0]) %>% ungroup()
    g <- ggplot(d, aes(AVISITN, pct, colour = COHORT, text = paste0(COHORT, "<br>Week ", AVISITN, ": n = ", n))) +
      geom_line(aes(group = COHORT)) + geom_point() + scale_colour_manual(values = pal) +
      scale_x_continuous(breaks = weeks()) + labs(x = "Week", y = "% of baseline participants", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  output$t_baseline <- renderDT({
    dat() %>% filter(AVISITN == 0) %>% group_by(Cohort = COHORT) %>%
      summarise(N = n(), `Age` = sprintf("%.1f (%.1f)", mean(AGEYR), sd(AGEYR)),
                `Education` = sprintf("%.1f (%.1f)", mean(EDCCNTU, na.rm = TRUE), sd(EDCCNTU, na.rm = TRUE)),
                `Female %` = sprintf("%.0f", 100 * mean(SEX == "Female", na.rm = TRUE)),
                `APOE4 carrier %` = sprintf("%.0f", 100 * mean(APOE4 == "carrier", na.rm = TRUE)),
                `Centiloids` = if (all(is.na(AMYLCENT))) "NA" else sprintf("%.1f (%.1f)", mean(AMYLCENT, na.rm = TRUE), sd(AMYLCENT, na.rm = TRUE)),
                `Baseline PACC` = sprintf("%.2f (%.2f)", mean(BASE, na.rm = TRUE), sd(BASE, na.rm = TRUE))) %>%
      datatable(options = list(dom = "t"), rownames = FALSE, caption = "Mean (SD) unless noted")
  })

  # --- trajectories
  output$p_spaghetti <- renderPlotly({
    d <- dat() %>% filter(!is.na(.data[[input$strat]]))
    req(nrow(d) > 0)
    set.seed(7); samp <- sample(unique(d$BID), min(120, n_distinct(d$BID)))
    m <- d %>% group_by(grp = .data[[input$strat]], AVISITN) %>% summarise(m = mean(AVAL), .groups = "drop")
    g <- ggplot() +
      geom_line(data = filter(d, BID %in% samp), aes(AVISITN, AVAL, group = BID), colour = "grey75", alpha = 0.5) +
      geom_line(data = m, aes(AVISITN, m, colour = grp), linewidth = 1.2) +
      geom_point(data = m, aes(AVISITN, m, colour = grp, text = sprintf("%s<br>Week %g: mean PACC %.2f", grp, AVISITN, m))) +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = weeks()) +
      labs(x = "Week", y = "PACC composite", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  output$p_meanchg <- renderPlotly({
    d <- dat() %>% filter(AVISITN > 0, !is.na(CHG), !is.na(.data[[input$strat]])) %>%
      group_by(grp = .data[[input$strat]], AVISITN) %>% summarise(m = mean(CHG), se = sd(CHG) / sqrt(n()), n = n(), .groups = "drop")
    req(nrow(d) > 0)
    g <- ggplot(d, aes(AVISITN, m, colour = grp, text = sprintf("%s<br>Week %g: %.2f +/- %.2f (n=%d)", grp, AVISITN, m, se, n))) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 5) + geom_line(aes(group = grp)) + geom_point() +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = weeks()) +
      labs(x = "Week", y = "Mean change from baseline", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  # --- MMRM
  observeEvent(input$fit, {
    d <- complete_cases()
    dropped <- setdiff(input$cohorts, levels(d$COHORT))
    if (nlevels(d$COHORT) < 2) { showNotification("At least two cohorts with post-baseline data are required.", type = "error"); return() }
    if (!isTRUE(input$ref %in% levels(d$COHORT))) {
      showNotification("The reference cohort has no analyzable post-baseline data in the current selection.", type = "error"); return() }
    d$COHORT <- relevel(d$COHORT, input$ref)
    res <- tryCatch(withProgress(message = "Fitting MMRM (unstructured covariance)...", {
      f <- mmrm(CHG ~ COHORT * AVISIT + BASE + AGEYR + EDCCNTU + APOE4 + us(AVISIT | BID), data = d,
                method = "Kenward-Roger", vcov = "Kenward-Roger-Linear")
      em <- emmeans(f, ~ COHORT | AVISIT)
      list(fit = f, n = n_distinct(d$BID), rows = nrow(d), dropped = dropped,
           lsm = as.data.frame(em) %>% left_join(distinct(d, AVISIT, week = AVISITN), by = "AVISIT"),
           diff = as.data.frame(contrast(em, "trt.vs.ctrl", adjust = "none")) %>%
             mutate(lower = estimate - qt(0.975, df) * SE, upper = estimate + qt(0.975, df) * SE))
    }), error = function(e) { showNotification(paste("MMRM fit failed:", conditionMessage(e)), type = "error", duration = 10); NULL })
    fit_res(res)
  })
  output$fit_msg <- renderPrint({
    f <- req(fit_res())
    cat(sprintf("Subjects: %d\nPost-baseline records: %d\nConverged: %s\nREML AIC: %.1f\n",
                f$n, f$rows, if (isTRUE(f$fit$opt_details$convergence == 0)) "yes" else "check", AIC(f$fit)))
    if (length(f$dropped)) cat("Excluded (no analyzable post-baseline data):", paste(f$dropped, collapse = ", "), "\n")
  })
  output$p_lsm <- renderPlotly({
    l <- req(fit_res())$lsm
    g <- ggplot(l, aes(week, emmean, colour = COHORT, text = sprintf("%s<br>Week %g: %.2f (%.2f, %.2f)", COHORT, week, emmean, lower.CL, upper.CL))) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 5) + geom_line(aes(group = COHORT)) + geom_point(size = 2) +
      scale_colour_manual(values = pal) + scale_x_continuous(breaks = weeks()) +
      labs(x = "Week", y = "LS-mean change from baseline", colour = NULL) + theme_minimal()
    ggplotly(g, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.25))
  })
  output$t_diff <- renderDT({
    req(fit_res())$diff %>% transmute(Contrast = contrast, Visit = AVISIT, Estimate = round(estimate, 3), SE = round(SE, 3),
                             df = round(df, 1), `95% CI` = sprintf("(%.3f, %.3f)", lower, upper), `p (unadjusted)` = signif(p.value, 2)) %>%
      datatable(options = list(pageLength = 12, dom = "tp"), rownames = FALSE,
                caption = "Two-sided t-tests with Kenward-Roger df; no multiplicity adjustment across visits or contrasts.")
  })

  # --- covariance comparison
  observeEvent(input$fitcov, {
    d <- complete_cases()
    if (nlevels(d$COHORT) < 2) { showNotification("At least two cohorts with post-baseline data are required.", type = "error"); return() }
    if (nrow(d) > 20000) { showNotification("Covariance comparison is limited to 20,000 post-baseline records; narrow the selection.", type = "error"); return() }
    mm <- CHG ~ COHORT * AVISIT + BASE + AGEYR + EDCCNTU + APOE4
    tv <- ~ as.integer(AVISIT) | BID
    res <- tryCatch(withProgress(message = "Fitting four covariance structures...", {
      fits <- list(
        "Compound symmetry"    = gls(mm, d, correlation = corCompSymm(form = tv), method = "REML"),
        "AR(1)"                = gls(mm, d, correlation = corAR1(form = tv), method = "REML"),
        "Heterogeneous AR(1)"  = gls(mm, d, correlation = corAR1(form = tv), weights = varIdent(form = ~ 1 | AVISIT), method = "REML"),
        "Unstructured"         = gls(mm, d, correlation = corSymm(form = tv), weights = varIdent(form = ~ 1 | AVISIT),
                                     method = "REML", control = glsControl(maxIter = 200, msMaxIter = 200)))
      tab <- tibble(Structure = names(fits), `Cov. parameters` = sapply(fits, function(f) attr(logLik(f), "df") - length(coef(f))),
                    AIC = round(sapply(fits, AIC), 1), BIC = round(sapply(fits, BIC), 1), logLik = round(sapply(fits, logLik), 1)) %>% arrange(AIC)
      list(tab = tab, lrt = anova(fits[["Heterogeneous AR(1)"]], fits[["Unstructured"]]))
    }), error = function(e) { showNotification(paste("Covariance fit failed:", conditionMessage(e)), type = "error", duration = 10); NULL })
    cov_res(res)
  })
  output$t_cov <- renderDT(datatable(req(cov_res())$tab, options = list(dom = "t"), rownames = FALSE))
  output$cov_lrt <- renderPrint({ a <- req(cov_res())$lrt
    cat(sprintf("LRT heterogeneous AR(1) vs unstructured: chi-sq = %.1f on %d df, p = %.3g\n",
                a$L.Ratio[2], a$df[2] - a$df[1], a$`p-value`[2])) })

  # --- validation
  output$t_recon <- renderDT({
    validate(need(!is.null(RECON), "Reconciliation table not bundled with this deployment."))
    RECON %>% transmute(Visit = AVISIT,
                        `R estimate` = signif(r_est, 6), `SAS estimate` = signif(sas_est, 6),
                        `R SE` = signif(r_se, 5), `SAS SE` = signif(sas_se, 5),
                        `R df` = round(r_df, 2), `SAS df` = round(sas_df, 2),
                        `|diff| estimate` = formatC(d_est, format = "e", digits = 2),
                        `|diff| SE` = formatC(d_se, format = "e", digits = 2),
                        `|diff| df` = round(d_df, 3)) %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })
  output$recon_summary <- renderPrint({
    req(RECON)
    cat(sprintf("Maximum absolute difference across %d visits: estimate = %.2e, SE = %.2e, df = %.2f\n",
                nrow(RECON), max(RECON$d_est), max(RECON$d_se), max(RECON$d_df)))
  })

  output$t_data <- renderDT(datatable(dat() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
                                      options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE, filter = "top"))
}

shinyApp(ui, server)

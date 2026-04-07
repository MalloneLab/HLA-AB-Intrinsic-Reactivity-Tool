### author: "Anna Mallone"
### date: "2025-12-10"

###############################################
## HLA-Antibodies IRT
## HLA Antibodies Intrinsic Reactivity Tool
###############################################

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(readr)
library(DT)
library(stringr)

###############################################
## Helper UI functions
###############################################

workflow_icon <- function(filename, fallback_icon, height = "70px") {
  if (file.exists(file.path("www", filename))) {
    tags$img(src = filename, height = height)
  } else {
    div(class = "icon-holder", icon(fallback_icon))
  }
}

home_icon_link <- function() {
  tags$a(
    href = "#",
    onclick = "Shiny.setInputValue('nav_target', 'home', {priority: 'event'})",
    style = "text-decoration:none;",
    if (file.exists(file.path("www", "home.png"))) {
      tags$img(
        src = "home.png",
        height = "60px",
        style = "cursor:pointer;"
      )
    } else {
      icon("home", style = "font-size: 50px; color: #0d6efd;")
    }
  )
}

page_header <- function(title_text) {
  fluidRow(
    column(
      width = 2,
      div(
        style = "padding-top: 5px;",
        home_icon_link()
      )
    ),
    column(
      width = 7,
      h2(title_text, style = "margin-top: 10px;")
    ),
    column(
      width = 3,
      div(
        style = "text-align: right;",
        if (file.exists(file.path("www", "logo.png"))) {
          tags$img(src = "logo.png", height = "70px")
        } else {
          NULL
        }
      )
    )
  )
}

###############################################
## Load and deduplicate reference table
###############################################

load_unspec_table <- function(path = "unspecificity_table.csv") {
  if (!file.exists(path)) {
    stop("Reference file 'unspecificity_table.csv' was not found in the application directory.")
  }
  
  raw_unspec_tbl <- readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      antigen = as.character(antigen),
      Class = as.factor(Class),
      unspec_prob = as.numeric(unspec_prob)
    )
  
  required_cols <- c(
    "antigen", "Class", "mean_freq_nonImm", "discordance_Imm",
    "disappearance_rate", "n_lots", "unspec_prob"
  )
  
  missing_cols <- setdiff(required_cols, names(raw_unspec_tbl))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Reference file is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  raw_unspec_tbl %>%
    group_by(antigen, Class) %>%
    summarise(
      mean_freq_nonImm   = mean(mean_freq_nonImm, na.rm = TRUE),
      discordance_Imm    = mean(discordance_Imm, na.rm = TRUE),
      disappearance_rate = mean(disappearance_rate, na.rm = TRUE),
      n_lots             = suppressWarnings(max(n_lots, na.rm = TRUE)),
      unspec_prob        = mean(unspec_prob, na.rm = TRUE),
      .groups = "drop"
    )
}

unspec_tbl <- tryCatch(
  load_unspec_table(),
  error = function(e) {
    data.frame(
      antigen = character(),
      Class = factor(),
      mean_freq_nonImm = numeric(),
      discordance_Imm = numeric(),
      disappearance_rate = numeric(),
      n_lots = numeric(),
      unspec_prob = numeric()
    )
  }
)

classI_antigens <- unspec_tbl %>%
  filter(Class == "I") %>%
  pull(antigen) %>%
  sort()

classII_antigens <- unspec_tbl %>%
  filter(Class == "II") %>%
  pull(antigen) %>%
  sort()

###############################################
## Helper functions
###############################################

safe_round <- function(x, digits = 3) {
  ifelse(is.na(x), NA_real_, round(x, digits))
}

detect_hla_class <- function(antigen) {
  ifelse(grepl("^(DR|DQ|DP)", antigen), "II", "I")
}

classify_unspec <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "no_data",
    p >= 0.5 ~ "high",
    p >= 0.3 ~ "intermediate",
    TRUE ~ "low"
  )
}

reco_text <- function(p) {
  if (is.na(p)) {
    return("No model information available for this antigen.")
  }
  if (p >= 0.5) {
    return("HIGH intrinsic reactivity probability (≥ 0.5). Orthogonal confirmation (e.g. Immucor Lifecodes) is strongly recommended before assigning this antigen as unacceptable.")
  }
  if (p >= 0.3) {
    return("INTERMEDIATE intrinsic reactivity probability (0.3–0.5). Interpret One Lambda positivity with caution; consider orthogonal confirmation when donor-specific or clinically relevant.")
  }
  return("LOW intrinsic reactivity probability (< 0.3). The One Lambda signal is generally more reliable, but should always be interpreted in the appropriate clinical context.")
}

category_fill <- function(cat) {
  dplyr::case_when(
    cat == "high" ~ "#f8d7da",
    cat == "intermediate" ~ "#fff3cd",
    cat == "low" ~ "#d4edda",
    TRUE ~ "#e2e3e5"
  )
}

percentile_text <- function(p, ref_probs) {
  if (is.na(p)) return("No reference match available")
  
  ref_probs <- ref_probs[!is.na(ref_probs)]
  if (length(ref_probs) == 0) return("No reference distribution available")
  
  top_pct <- max(1, ceiling(mean(ref_probs >= p) * 100))
  paste0("This antigen is in the top ", top_pct, "% highest intrinsic reactivity antigens")
}

extract_hla_antigens <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(character(0))
  
  txt_clean <- gsub("[\r\n\t]", " ", txt)
  
  hits <- unlist(regmatches(
    txt_clean,
    gregexpr(
      "\\b(?:A\\d{1,2}|B\\d{1,2}|Cw\\d{1,2}|DR\\d{1,2}|DQ\\d{1,2}|DP\\d{1,2})\\b",
      txt_clean,
      perl = TRUE
    )
  ))
  
  hits <- unique(trimws(hits))
  hits[hits != ""]
}

fusion_text_to_table <- function(txt, ref_tbl, patient_id = "pasted_case") {
  ag <- extract_hla_antigens(txt)
  
  if (length(ag) == 0) {
    return(
      data.frame(
        patient_id = character(),
        antigen = character(),
        Class = character(),
        intrinsic_reactivity_probability = numeric(),
        percentile_interpretation = character(),
        category = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  ref_probs <- ref_tbl$unspec_prob
  
  data.frame(
    patient_id = patient_id,
    antigen = ag,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Class = detect_hla_class(antigen)
    ) %>%
    left_join(
      ref_tbl %>% mutate(Class = as.character(Class)),
      by = c("antigen", "Class")
    ) %>%
    mutate(
      intrinsic_reactivity_probability = safe_round(unspec_prob),
      percentile_interpretation = vapply(unspec_prob, percentile_text, character(1), ref_probs = ref_probs),
      category = classify_unspec(unspec_prob)
    ) %>%
    select(
      patient_id,
      antigen,
      Class,
      intrinsic_reactivity_probability,
      percentile_interpretation,
      category
    ) %>%
    arrange(desc(ifelse(is.na(intrinsic_reactivity_probability), -1, intrinsic_reactivity_probability)))
}

traffic_counts <- function(df, prob_col = "intrinsic_reactivity_probability") {
  if (nrow(df) == 0) {
    return(data.frame(category = character(), n = integer()))
  }
  
  p <- df[[prob_col]]
  data.frame(category = classify_unspec(p)) %>%
    count(category, name = "n") %>%
    right_join(
      data.frame(category = c("high", "intermediate", "low", "no_data"), stringsAsFactors = FALSE),
      by = "category"
    ) %>%
    mutate(n = ifelse(is.na(n), 0L, n))
}

###############################################
## Theme
###############################################

app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Roboto"),
  heading_font = font_google("Roboto Slab")
)

###############################################
## UI
###############################################

ui <- fluidPage(
  theme = app_theme,
  
  tags$head(
    tags$style(HTML("
      .small-note { font-size: 0.88rem; color: #555; }
      .box-panel {
        background-color: #ffffff;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 18px;
        margin-bottom: 18px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.05);
      }
      .info-box {
        background-color: #f8f9fa;
        border-left: 4px solid #d6d9dc;
        padding: 10px 12px;
        margin-bottom: 16px;
        border-radius: 4px;
        font-size: 0.85rem;
        color: #555;
      }
      .traffic-box {
        display: inline-block;
        min-width: 120px;
        padding: 12px;
        margin-right: 10px;
        margin-bottom: 10px;
        border-radius: 8px;
        text-align: center;
        font-weight: 600;
      }
      .traffic-title {
        font-size: 0.9rem;
        display: block;
      }
      .traffic-count {
        font-size: 1.4rem;
        display: block;
      }
      .home-card {
        background-color: white;
        border: 1px solid #dee2e6;
        border-radius: 12px;
        padding: 25px;
        text-align: center;
        cursor: pointer;
        transition: all 0.2s ease;
        height: 280px;
        overflow: hidden;
      }
      .home-card:hover {
        box-shadow: 0 6px 18px rgba(0,0,0,0.15);
        transform: translateY(-5px);
      }
      .home-card h3 {
        margin-top: 15px;
        font-size: 1.25rem;
      }
      .home-card p {
        font-size: 14px;
        color: #555;
        white-space: normal;
        overflow-wrap: break-word;
      }
      .icon-holder {
        font-size: 44px;
        margin-top: 5px;
        color: #0d6efd;
      }
      .center-subtitle {
        text-align: center;
        font-size: 1rem;
        color: #666;
        margin-top: -8px;
      }
      .center-description {
        text-align: center;
        font-size: 1rem;
        margin-top: 20px;
        margin-bottom: 20px;
      }
      img {
        transition: transform 0.2s ease, opacity 0.2s ease;
      }
      img:hover {
        transform: scale(1.04);
        opacity: 0.95;
      }
    "))
  ),
  
  tabsetPanel(
    id = "mainnav",
    type = "hidden",
    
    ###############################################
    ## HOME
    ###############################################
    tabPanel(
      title = "home",
      value = "home",
      
      br(),
      
      div(
        class = "info-box",
        "This tool is intended for research and educational use only. It does not replace expert laboratory interpretation or clinical decision-making."
      ),
      
      if (file.exists(file.path("www", "logo_extended.png"))) {
        div(
          style = "text-align:center;",
          tags$img(src = "logo_extended.png", height = "90px")
        )
      } else {
        tagList(
          h1("HLA-Antibodies IRT", align = "center"),
          div(class = "center-subtitle", "HLA Antibodies Intrinsic Reactivity Tool")
        )
      },
      
      div(
        class = "center-description",
        "This tool enables exploration and interpretation of assay-intrinsic reactivity in HLA antibody testing. Select one of the workflows below."
      ),
      
      br(),
      
      fluidRow(
        column(
          width = 4,
          div(
            class = "home-card",
            onclick = "Shiny.setInputValue('nav_target', 'antigen', {priority: 'event'})",
            workflow_icon("lookup.png", "search", height = "85px"),
            h3("Antigen Lookup"),
            p("Query one antigen and inspect its intrinsic reactivity profile.")
          )
        ),
        column(
          width = 4,
          div(
            class = "home-card",
            onclick = "Shiny.setInputValue('nav_target', 'singlefusion', {priority: 'event'})",
            workflow_icon("pastefusion.png", "clipboard", height = "85px"),
            h3("Query Single Fusion Assignment"),
            p("Copy-paste Fusion final assignments and obtain color-coded intrinsic reactivity profile.")
          )
        ),
        column(
          width = 4,
          div(
            class = "home-card",
            onclick = "Shiny.setInputValue('nav_target', 'multipatient', {priority: 'event'})",
            workflow_icon("patients.png", "users", height = "85px"),
            h3("Multiple Patient Reports"),
            p("Query multiple patients assignments and generate combined reports.")
          )
        )
      )
    ),
    
    ###############################################
    ## ANTIGEN LOOKUP
    ###############################################
    tabPanel(
      title = "antigen",
      value = "antigen",
      
      br(),
      page_header("Antigen Lookup"),
      
      div(
        class = "box-panel",
        selectInput("antigen_I", "Class I antigen",
                    choices = c("None" = "", classI_antigens)),
        selectInput("antigen_II", "Class II antigen",
                    choices = c("None" = "", classII_antigens)),
        p(
          class = "small-note",
          "Select one antigen. If both Class I and Class II are selected, the Class II antigen will be prioritised."
        ),
        hr(),
        div(
          class = "small-note",
          "This section summarizes intrinsic reactivity metrics based on:",
          tags$ul(
            tags$li("Background reactivity in non-immunized male sera"),
            tags$li("Cross-platform One Lambda/Immucor discordance"),
            tags$li("Post-transplant disappearance dynamics")
          )
        ),
        br(),
        h4("Antigen-level intrinsic reactivity profile"),
        div(
          strong("Columns: "),
          "antigen ", actionLink("info_a", "(i)"), " · ",
          "Class ", actionLink("info_c", "(i)"), " · ",
          "mean_freq_nonImm ", actionLink("info_mn", "(i)"), " · ",
          "discordance_Imm ", actionLink("info_di", "(i)"), " · ",
          "disappearance_rate ", actionLink("info_dr", "(i)"), " · ",
          "n_lots ", actionLink("info_nl", "(i)"), " · ",
          "unspec_prob ", actionLink("info_up", "(i)")
        ),
        br(), br(),
        DTOutput("antigen_table"),
        br(),
        h4("Interpretation"),
        uiOutput("antigen_reco_ui"),
        br(),
        h4("Traffic-light summary"),
        uiOutput("antigen_traffic"),
        br(),
        h4("Relative position"),
        uiOutput("antigen_percentile_ui")
      )
    ),
    
    ###############################################
    ## SINGLE FUSION QUERY
    ###############################################
    tabPanel(
      title = "singlefusion",
      value = "singlefusion",
      
      br(),
      page_header("Query Single Fusion Assignment"),
      
      div(
        class = "box-panel",
        textInput("single_patient_id", "Patient/sample ID", value = "single_case_1"),
        textAreaInput(
          "single_fusion_paste",
          label = NULL,
          placeholder = "Paste here the text copied from Fusion via 'Copy final assignment'...",
          rows = 8,
          width = "100%"
        ),
        fluidRow(
          column(6, actionButton("calculate_single", "Calculate")),
          column(6, actionButton("clear_single", "Clear assignments"))
        ),
        br(),
        p(
          class = "small-note",
          "The app will automatically extract recognizable HLA antigen assignments ",
          "(e.g. A80, B76, Cw14, DR53, DQ7, DP1), match them to the intrinsic reactivity table, ",
          "and display a color-coded interpretation."
        ),
        br(),
        downloadButton("download_single_report", "Download report"),
        br(), br(),
        h4("Traffic-light summary"),
        uiOutput("single_fusion_traffic"),
        br(),
        h4("Parsed assignments"),
        DTOutput("single_fusion_table")
      )
    ),
    
    ###############################################
    ## MULTIPLE PATIENT REPORTS
    ###############################################
    tabPanel(
      title = "multipatient",
      value = "multipatient",
      
      br(),
      page_header("Multiple Patient Reports"),
      
      div(
        class = "box-panel",
        textInput("multi_patient_id", "Patient/sample ID", value = "patient_1"),
        textAreaInput(
          "multi_fusion_paste",
          label = "Paste Fusion final assignment",
          placeholder = "Paste here the text copied from Fusion via 'Copy final assignment'...",
          rows = 8,
          width = "100%"
        ),
        fluidRow(
          column(4, actionButton("add_multi_patient", "Add patient")),
          column(4, actionButton("clear_multi_input", "Clear current text")),
          column(4, actionButton("clear_all_multi", "Clear all patients"))
        ),
        br(),
        p(
          class = "small-note",
          "Add patients one by one by entering a patient ID and pasting the corresponding Fusion final assignment."
        ),
        br(),
        downloadButton("download_multi_report", "Download combined report"),
        br(), br(),
        h4("Traffic-light overview by patient"),
        uiOutput("multi_patient_summary_ui"),
        br(),
        h4("Selected patient report"),
        DTOutput("multi_patient_table")
      )
    )
  )
)

###############################################
## SERVER
###############################################

server <- function(input, output, session) {
  
  observeEvent(input$nav_target, {
    updateTabsetPanel(session, "mainnav", selected = input$nav_target)
  })
  
  ###############################################
  ## Antigen selection
  ###############################################
  
  selected_antigen <- reactive({
    if (!is.null(input$antigen_II) && nzchar(input$antigen_II)) return(input$antigen_II)
    if (!is.null(input$antigen_I)  && nzchar(input$antigen_I))  return(input$antigen_I)
    NULL
  })
  
  antigen_data <- reactive({
    req(selected_antigen())
    unspec_tbl %>% filter(antigen == selected_antigen())
  })
  
  ###############################################
  ## Info popups
  ###############################################
  
  observeEvent(input$info_a, {
    showModal(modalDialog("HLA antigen specificity (for example DR53, DP1, or A80).", easyClose = TRUE))
  })
  observeEvent(input$info_c, {
    showModal(modalDialog("HLA Class I (A/B/C) or Class II (DR/DQ/DP).", easyClose = TRUE))
  })
  observeEvent(input$info_mn, {
    showModal(modalDialog("Mean One Lambda positivity frequency in the non-immunized male reference cohort across available lots (range 0–1).", easyClose = TRUE))
  })
  observeEvent(input$info_di, {
    showModal(modalDialog("Proportion of baseline One Lambda-positive signals not confirmed by Immucor (OL+/Imm−).", easyClose = TRUE))
  })
  observeEvent(input$info_dr, {
    showModal(modalDialog("Proportion of OL+/Imm− signals that fall below the One Lambda cutoff by 12 months post-transplant.", easyClose = TRUE))
  })
  observeEvent(input$info_nl, {
    showModal(modalDialog("Number of One Lambda lots available for this antigen.", easyClose = TRUE))
  })
  observeEvent(input$info_up, {
    showModal(modalDialog("Model-predicted probability that a One Lambda-positive signal reflects assay-intrinsic reactivity.", easyClose = TRUE))
  })
  
  ###############################################
  ## Antigen lookup outputs
  ###############################################
  
  output$antigen_table <- renderDT({
    req(nrow(antigen_data()) > 0)
    
    df <- antigen_data() %>%
      mutate(
        mean_freq_nonImm = safe_round(mean_freq_nonImm),
        discordance_Imm = safe_round(discordance_Imm),
        disappearance_rate = safe_round(disappearance_rate),
        unspec_prob = safe_round(unspec_prob)
      )
    
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
  
  output$antigen_reco_ui <- renderUI({
    req(nrow(antigen_data()) > 0)
    
    df <- antigen_data()
    p <- df$unspec_prob[1]
    
    div(
      style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 12px;",
      tags$p(tags$strong("Antigen: "), paste0(df$antigen[1], " (Class ", df$Class[1], ")")),
      tags$p(tags$strong("Intrinsic reactivity probability: "), sprintf("%.3f", p)),
      tags$p(reco_text(p))
    )
  })
  
  output$antigen_traffic <- renderUI({
    req(nrow(antigen_data()) > 0)
    
    df <- antigen_data()
    p <- df$unspec_prob[1]
    cat <- classify_unspec(p)
    
    tags$div(
      class = "traffic-box",
      style = paste0("background:", category_fill(cat), ";"),
      tags$span(class = "traffic-title", toupper(cat)),
      tags$span(class = "traffic-count", sprintf("%.3f", p))
    )
  })
  
  output$antigen_percentile_ui <- renderUI({
    req(nrow(antigen_data()) > 0)
    
    df <- antigen_data()
    p <- df$unspec_prob[1]
    
    div(
      style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 12px;",
      percentile_text(p, unspec_tbl$unspec_prob)
    )
  })
  
  ###############################################
  ## Single Fusion assignment
  ###############################################
  
  single_store <- reactiveVal(
    data.frame(
      patient_id = character(),
      antigen = character(),
      Class = character(),
      intrinsic_reactivity_probability = numeric(),
      percentile_interpretation = character(),
      category = character(),
      stringsAsFactors = FALSE
    )
  )
  
  observeEvent(input$calculate_single, {
    validate(
      need(nzchar(trimws(input$single_fusion_paste)), "Please paste Fusion final assignment text.")
    )
    
    pid <- trimws(input$single_patient_id)
    if (!nzchar(pid)) pid <- "single_case_1"
    
    parsed <- fusion_text_to_table(
      txt = input$single_fusion_paste,
      ref_tbl = unspec_tbl,
      patient_id = pid
    )
    
    single_store(parsed)
  })
  
  observeEvent(input$clear_single, {
    single_store(
      data.frame(
        patient_id = character(),
        antigen = character(),
        Class = character(),
        intrinsic_reactivity_probability = numeric(),
        percentile_interpretation = character(),
        category = character(),
        stringsAsFactors = FALSE
      )
    )
    updateTextAreaInput(session, "single_fusion_paste", value = "")
  })
  
  output$single_fusion_traffic <- renderUI({
    cnt <- traffic_counts(single_store(), "intrinsic_reactivity_probability")
    
    tags$div(
      lapply(seq_len(nrow(cnt)), function(i) {
        tags$div(
          class = "traffic-box",
          style = paste0("background:", category_fill(cnt$category[i]), ";"),
          tags$span(class = "traffic-title", toupper(cnt$category[i])),
          tags$span(class = "traffic-count", cnt$n[i])
        )
      })
    )
  })
  
  output$single_fusion_table <- renderDT({
    df <- single_store() %>%
      select(
        antigen,
        Class,
        percentile_interpretation,
        category
      ) %>%
      rename(
        Antigen = antigen,
        Class = Class,
        `Percentile interpretation` = percentile_interpretation,
        Category = category
      )
    
    datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        columnDefs = list(
          list(width = '45%', targets = 2)
        )
      )
    ) %>%
      formatStyle(
        "Category",
        backgroundColor = styleEqual(
          c("high", "intermediate", "low", "no_data"),
          c("#f8d7da", "#fff3cd", "#d4edda", "#e2e3e5")
        )
      ) %>%
      formatStyle(
        "Percentile interpretation",
        `white-space` = "normal"
      )
  })
  
  output$download_single_report <- downloadHandler(
    filename = function() paste0("single_fusion_intrinsic_reactivity_report_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- single_store()
      
      summary_tbl <- df %>%
        count(category, name = "n_antigens")
      
      report_tbl <- df %>%
        select(patient_id, antigen, Class, percentile_interpretation, category)
      
      wb <- createWorkbook()
      
      addWorksheet(wb, "Traffic-light summary")
      writeData(wb, "Traffic-light summary", summary_tbl)
      
      addWorksheet(wb, "Parsed assignments")
      writeData(wb, "Parsed assignments", report_tbl)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  ###############################################
  ## Multiple Patient Reports
  ###############################################
  
  multi_store <- reactiveVal(
    data.frame(
      patient_id = character(),
      antigen = character(),
      Class = character(),
      intrinsic_reactivity_probability = numeric(),
      percentile_interpretation = character(),
      category = character(),
      stringsAsFactors = FALSE
    )
  )
  
  observeEvent(input$add_multi_patient, {
    validate(
      need(nzchar(trimws(input$multi_fusion_paste)), "Please paste Fusion final assignment text.")
    )
    
    pid <- trimws(input$multi_patient_id)
    if (!nzchar(pid)) pid <- paste0("patient_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    
    parsed <- fusion_text_to_table(
      txt = input$multi_fusion_paste,
      ref_tbl = unspec_tbl,
      patient_id = pid
    )
    
    multi_store(
      multi_store() %>%
        bind_rows(parsed) %>%
        distinct()
    )
  })
  
  observeEvent(input$clear_multi_input, {
    updateTextAreaInput(session, "multi_fusion_paste", value = "")
  })
  
  observeEvent(input$clear_all_multi, {
    multi_store(
      data.frame(
        patient_id = character(),
        antigen = character(),
        Class = character(),
        intrinsic_reactivity_probability = numeric(),
        percentile_interpretation = character(),
        category = character(),
        stringsAsFactors = FALSE
      )
    )
    updateTextAreaInput(session, "multi_fusion_paste", value = "")
  })
  
  multi_summary_tbl <- reactive({
    df <- multi_store()
    
    if (nrow(df) == 0) {
      return(
        data.frame(
          patient_id = character(),
          n_high = integer(),
          n_intermediate = integer(),
          n_low = integer(),
          n_no_data = integer(),
          stringsAsFactors = FALSE
        )
      )
    }
    
    df %>%
      group_by(patient_id) %>%
      summarise(
        n_high = sum(category == "high", na.rm = TRUE),
        n_intermediate = sum(category == "intermediate", na.rm = TRUE),
        n_low = sum(category == "low", na.rm = TRUE),
        n_no_data = sum(category == "no_data", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(patient_id)
  })
  
  selected_multi_patient <- reactiveVal(NULL)
  
  observeEvent(input$select_multi_patient, {
    selected_multi_patient(input$select_multi_patient)
  })
  
  output$multi_patient_summary_ui <- renderUI({
    df <- multi_summary_tbl()
    
    if (nrow(df) == 0) {
      return(tags$p(class = "small-note", "No patients added yet."))
    }
    
    tagList(
      lapply(seq_len(nrow(df)), function(i) {
        pid <- df$patient_id[i]
        
        tags$div(
          style = "margin-bottom: 14px; padding: 10px 12px; border: 1px solid #dee2e6; border-radius: 8px; background-color: #ffffff;",
          
          tags$div(
            style = "margin-bottom: 8px;",
            tags$a(
              href = "#",
              onclick = sprintf("Shiny.setInputValue('select_multi_patient', '%s', {priority: 'event'})", pid),
              style = "font-weight: 700; font-size: 16px; text-decoration: none;",
              pid
            )
          ),
          
          tags$div(
            style = "display:flex; gap:10px; flex-wrap:wrap;",
            
            tags$div(
              class = "traffic-box",
              style = paste0("background:", category_fill("low"), ";"),
              tags$span(class = "traffic-title", "LOW"),
              tags$span(class = "traffic-count", df$n_low[i])
            ),
            
            tags$div(
              class = "traffic-box",
              style = paste0("background:", category_fill("high"), ";"),
              tags$span(class = "traffic-title", "HIGH"),
              tags$span(class = "traffic-count", df$n_high[i])
            ),
            
            tags$div(
              class = "traffic-box",
              style = paste0("background:", category_fill("intermediate"), ";"),
              tags$span(class = "traffic-title", "INTERMEDIATE"),
              tags$span(class = "traffic-count", df$n_intermediate[i])
            ),
            
            tags$div(
              class = "traffic-box",
              style = paste0("background:", category_fill("no_data"), ";"),
              tags$span(class = "traffic-title", "NO DATA"),
              tags$span(class = "traffic-count", df$n_no_data[i])
            )
          )
        )
      })
    )
  })
  
  output$multi_patient_table <- renderDT({
    req(selected_multi_patient())
    
    detail <- multi_store() %>%
      filter(patient_id == selected_multi_patient()) %>%
      select(
        patient_id,
        antigen,
        Class,
        percentile_interpretation
      ) %>%
      rename(
        `Patient ID` = patient_id,
        `Antigen` = antigen,
        `Class` = Class,
        `Percentile interpretation` = percentile_interpretation
      )
    
    datatable(
      detail,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        autoWidth = TRUE,
        columnDefs = list(
          list(width = '45%', targets = 3)
        )
      )
    ) %>%
      formatStyle(
        "Percentile interpretation",
        `white-space` = "normal"
      )
  })
  
  output$download_multi_report <- downloadHandler(
    filename = function() paste0("multiple_patient_intrinsic_reactivity_report_", Sys.Date(), ".xlsx"),
    content = function(file) {
      summary_tbl <- multi_summary_tbl()
      
      detail_tbl <- multi_store() %>%
        select(
          patient_id,
          antigen,
          Class,
          percentile_interpretation
        ) %>%
        arrange(patient_id, antigen)
      
      wb <- createWorkbook()
      
      addWorksheet(wb, "Patient summary")
      writeData(wb, "Patient summary", summary_tbl)
      
      addWorksheet(wb, "Detailed report")
      writeData(wb, "Detailed report", detail_tbl)
      
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

###############################################
## RUN APP
###############################################

shinyApp(ui = ui, server = server)
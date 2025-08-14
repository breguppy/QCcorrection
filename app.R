library(shiny)
library(shinyBS)
library(bslib)
library(dplyr)
library(shinycssloaders)
library(openxlsx)
library(readxl)
library(tools)
library(ggplot2)
library(tidyr)
library(purrr)
library(tidyverse)
source("R/helpers.R")
source("R/processing_helpers.R")
source("R/met_scatter_rf.R")
source("R/met_scatter_loess.R")
source("R/download_helpers.R")
source("R/plotting_helpers.R")
source("R/plotting_rsd_comparisons.R")


ui <- fluidPage(
  theme = bs_theme(preset = "cosmo"),
  titlePanel("QC Correction for Metabolomics Data"),
  
  navset_tab(
    id = "main_steps",
    
    #--- Step 1: Import Raw Data
    nav_panel(
      title = "1. Import Raw Data",
      
      #-- Upload raw data
      card(full_screen = TRUE, # Allows user to open table into full screen
           layout_sidebar(
             sidebar = sidebar(
               #--- Upload raw data
               tags$h4("1.1 Upload Raw Data"),
               fileInput(
                 inputId = "file1",
                 label = "Choose Raw Data File (.csv, .xls, or .xlsx)",
                 accept = c(".csv", ".xls", ".xlsx"),
                 buttonLabel = "Browse...",
                 placeholder = "No file selected"
               ),
               tags$h6("Raw data must be on the first sheet of the .xls or .xlsx file."),
               tags$h6(
                 "When sorted by injection order, the data must begin and end with a QC samples."
               ),
               width = 400
             ),
             tags$div(style = "overflow-x: auto; overflow-y: auto; max-height: 400px;", tableOutput("contents"))
           )),
      #-- User provides basic information
      card(layout_sidebar(
        sidebar = sidebar(
          #--- select non-metabolite columns
          tags$h4("1.2 Select non-metabolite columns"),
          tags$h6("Please select columns for sample, batch, class, and order."),
          uiOutput("column_selectors"),
          uiOutput("column_warning"),
          tooltip(
            checkboxInput(
              inputId = "withhold_cols",
              label = "Withhold additional columns from correction?",
              value = FALSE
            ),
            "Check the box if there are more non-metabolite columns (other than what is listed above) or specific metabolite columns that should not be corrected.",
            placement = "right"
          ),
          uiOutput("n_withhold_ui"),
          uiOutput("withhold_selectors_ui"),
          width = 400,
        ),
        uiOutput("basic_info")
      ), ),
      #-- Filter data based on missing values
      card(layout_sidebar(
        sidebar = sidebar(
          #--- Filter metabolites
          tags$h4("1.3 Filter Raw Data"),
          tooltip(
            sliderInput(
              inputId = "Frule",
              label = "Acceptable % of missing values per metabolite",
              min = 0,
              max = 100,
              value = 20
            ),
            "Metabolites with more than the acceptable % of missing values will be removed from the data.",
            placement = "right"
          ),
          width = 400
        ),
        uiOutput("filter_info"),
      ), ),
      #-- move on to step 2 button
      card(
        actionButton(
          inputId = "next_correction",
          label = "Next: Choose Correction Settings",
          class = "btn-primary btn-lg"
        ),
      )
    ),
    
    #--- Step 2: Correction settings
    nav_panel(
      title = "2. Correction Settings",
      
      #-- User selects correction settings
      card(
        style = "background-color: #eeeeee;",
        tags$h4("2.1 Choose Correction Settings"),
        uiOutput("qc_missing_value_warning"),
        fluidRow(
          column(3, tags$h5("Impute Missing QC Values"), uiOutput("qcImpute")),
          column(
            3,
            tags$h5("Impute Missing Sample Values"),
            uiOutput("sampleImpute")
          ),
          column(
            3,
            #--- Choose Correction method
            tags$h5("Choose Correction Method"),
            uiOutput("correctionMethod"),
          ),
          column(
            3,
            tags$h5("Unavailable Options"),
            uiOutput("unavailable_options")
          ),
          actionButton(
            inputId = "correct",
            label = "Correct Data with Selected Settings",
            class = "btn-primary btn-lg"
          ),
        ),
      ),
      #-- User selects post-correction filtering and transformation/normalization
      card(
        layout_sidebar(
          sidebar = sidebar(
            # After correction filtering
            tags$h4("2.2 Post-Correction Filtering"),
            tooltip(
              checkboxInput(
                inputId = "remove_imputed",
                label = "Remove imputed values after correction?",
                value = FALSE
              ),
              "Check this box if you want to the corrected data to have the same missing values as the raw data.",
              placement = "right"
            ),
            tooltip(
              checkboxInput(
                inputId = "post_cor_filter",
                label = "Don't filter metabolites based on QC RSD%",
                value = FALSE
              ),
              "Check this box if you don't want any metabolites removed post-correction.",
              placement = "right"
            ),
            conditionalPanel(
              "input.post_cor_filter == false",
              tooltip(
                sliderInput(
                  inputId = "rsd_filter",
                  label = "Metabolite RSD% threshold for QC samples",
                  min = 0,
                  max = 100,
                  value = 50
                ),
                "Metabolites with QC RSD% above this value will be removed from the corrected data.",
                placement = "right"
              )
            ),
            tags$hr(),
            
            # After correction scaling / normalization
            tags$h4("2.3 Post-Correction Transformation or Normalization"),
            radioButtons(
              inputId = "transform",
              label = "Method",
              choices = list(
                "Log 2 Transformation" = "log2",
                "Total Ratiometically Normalized (TRN)" = "TRN",
                "None" = "none"
              ),
              selected = "none"
            ),
            tooltip(
              checkboxInput(
                inputId = "ex_ISTD",
                label = "Exclude Internal Standards from post-correction transformation/normalization.",
                value = TRUE
              ),
              "Check this box if you do not want internal standards to be included in the transformation or normalization calculation.",
              placement = "right"
            ),
            conditionalPanel(
              "input.transform == 'TRN'",
              tooltip(
                checkboxInput(
                  inputId = "trn_withhold_checkbox",
                  label = "Withold column(s) from TRN?",
                  value = FALSE
                ),
                "Check this box if there are any columns that should not count in TRN (i.e. TIC column). Sample, batch, class and order are already excluded.",
                placement = "right"
              )
            ),
            uiOutput("trn_withhold_ui"),
            uiOutput("trn_withhold_selectors_ui"),
            width = 400,
          ),
          tags$div(
            style = "overflow-x: auto; overflow-y: auto; max-height: 400px; border: 1px solid #ccc;",
            tableOutput("cor_data") %>% withSpinner(color = "#404040")
          ),
          uiOutput("post_cor_filter_info") %>% withSpinner(color = "#404040"),
        ),
      ),
      card(
        style = "background-color: #eeeeee;",
        fluidRow(column(
          4,
          tags$h4("2.3 Identify Control Group"),
          tooltip(
            checkboxInput(
              inputId = "no_control",
              label = "No control group.",
              value = FALSE
            ),
            "Check the box if The data does not have a control group.",
            placement = "right"
          ),
          conditionalPanel(
            "input.no_control == false",
            uiOutput("control_class_selector")
          )
        ),
        column(4,
               tags$h4("2.5 Download Corrected Data"),
               tooltip(
                 checkboxInput(
                   inputId = "keep_corrected_qcs",
                   label = "Include QCs in corrected data file.",
                   value = FALSE
                 ),
                 "Check the box if you want corrected QC values in the downloaded corrected data file.",
                 placement = "right"
               ),
               uiOutput("download_corr_btn", container = div, style = "position: absolute; bottom: 15px; right: 15px;")
        ))
      ),
      card(
        actionButton(
          inputId = "next_visualization",
          label = "Next: Evaluate and Visualize Correction",
          class = "btn-primary btn-lg"
        ),
      )
    ),
    
    #--- Step 3: plots
    nav_panel(
      title = "3. Evaluation Metrics and Visualization",
      
      card(layout_sidebar(
        # RSD Plots
        sidebar = sidebar(
          tags$h4("3.1 RSD Evaluation"),
          radioButtons(
            inputId = "rsd_cal",
            label = "Calculate RSD by",
            choices = list("Metabolite" = "met", "Class and Metabolite" = "class_met"),
            selected = "met"
          ),
          width = 400,
        ),
        plotOutput("rsd_comparison_plot", height = "500px", width = "1100px")
      )),
      card(layout_sidebar(
        sidebar = sidebar(
          #-- PCA plots
          tags$h4("3.2 PCA Evaluation"),
          radioButtons(
            inputId = "color_col",
            label = "Color PCA by",
            choices = list("batch" = "batch", "class" = "class"),
            selected = "batch"
          ),
          width = 400,
        ),
        plotOutput("pca_plot", height = "500px", width = "1100px")
      )),
      card(layout_sidebar(
        sidebar = sidebar(
          #--- plot metabolite
          tags$h4("3.4 Metabolite Scatter plots"),
          uiOutput("met_plot_selectors"),
          width = 400,
        ),
        plotOutput("metab_scatter", height = "600px", width = "600px"),
      )),
      card(layout_sidebar(
        #-- select figure format
        sidebar = sidebar(
          tooltip(
            radioButtons(
              inputId = "fig_format",
              label = "Select figure format:",
              choices = c("PDF" = "pdf", "PNG" = "png"),
              selected = "pdf"
            ),
            "All figures will be saved in this format after clicking download button.",
            placement = "right"
          ),
          width = 400,
        ),
        uiOutput(
          "download_fig_zip_btn"
          #container = div,
          #style = "position: absolute; bottom: 15px; right: 15px;"
        ),
        uiOutput("progress_ui"),
      )),
      card(
        actionButton(
          inputId = "next_export",
          label = "Next: Export Corrected Data and Plots",
          class = "btn-primary btn-lg"
        ),
      )
    ),
    #--- Step 5: Export Data
    nav_panel(title = "4. Export Corrected Data and Plots", card(
      card_title("Download Data and Plots"),
      tags$span("Download all to get a ", icon("folder"), " zipped folder containing:"),
      fluidRow(
        column(4, tags$span(icon("file-excel"), "corrected_data_today's date.xlsx"),
               tags$ul(
                 tags$li("0. Raw Data"),
                 tags$li("1. Correction Settings"),
                 tags$li("2. Drift Normalized"),
                 tags$li("3. Scaled or Normalized"),
                 tags$li("4. Grouped Data Organized"),
                 tags$li("5. Grouped Data Fold Change (only when provided a control class)"),
                 tags$li("Appendix1. Metaboanalyst Ready")
               )),
        column(4, tags$span(icon("folder"), " figures"),
               tags$ul(
                 tags$li(icon("folder"), " metabolite figures"),
                 tags$li(icon("folder"), " RSD figures"),
                 tags$li(icon("folder"), "PCA plots")
               ))
      ),
      uiOutput("download_all_ui")
    ))
  )
)

# Define server
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    message("Session ended; calling stopApp()")
    stopApp()
  })
  # Define reactive values so we can save a snapshot after each tab.
  # This ensure we have finalized versions for plotting after correction.
  rv <- reactiveValues(cleaned = NULL, 
                       filtered = NULL, 
                       imputed = NULL, 
                       filtered_corrected = NULL,
                       transformed = NULL,
                       params = NULL)
  
  #-- Import Raw Data
  data_raw <- reactive({
    req(input$file1)
    read_raw_data(input$file1$datapath)
  }) #%>% bindCache(input$file1$datapath)
  
  #–– View Raw Data
  output$contents <- renderTable(data_raw())
  
  #-- reactive for metadata column selection
  selections_r <- reactive({
    list(
      sample   = input$sample_col %||% "",
      batch    = input$batch_col  %||% "",
      class    = input$class_col  %||% "",
      order    = input$order_col  %||% ""
    )
  }) %>% debounce(200)
  #–– Column selection & warnings
  output$column_selectors <- renderUI({
    req(data_raw())
    cols <- names(data_raw())
    nonMetColSelectionUI(cols)
  })
  output$column_warning <- renderUI({
    req(data_raw())
    sel <- selections_r()
    selected <- c(sel$sample, sel$batch, sel$class, sel$order)
    columnWarningUI(data_raw(), selected)
  })
  
  #-- Option to withhold more columns from correction.
  withheld_ids_r <- reactive({
    if (!isTRUE(input$withhold_cols)) return(character(0))
    n <- input$n_withhold %||% 0
    if (n <= 0) return(character(0))
    paste0("withhold_col_", seq_len(n))
  })
  withheld_r <- reactive({
    ids <- withheld_ids_r()
    if (!length(ids)) return(character(0))
    vals <- vapply(ids, function(id) input[[id]] %||% "", character(1))
    vals <- unique(vals[nzchar(vals)])
    
    sel <- selections_r()
    metadata <- c(sel$sample, sel$batch, sel$class, sel$order)
    setdiff(vals, metadata)
  }) %>% debounce(200)
  observe({
    req(data_raw())
    
    max_withhold <- max(ncol(data_raw()) - 4, 0)
    
    output$n_withhold_ui <- renderUI({
      if (isTRUE(input$withhold_cols)) {
        numericInput(
          inputId = "n_withhold",
          label = "Number of columns to withold from correction",
          value = 1,
          min = 1,
          max = max_withhold
        )
      }
    })
  })
  output$withhold_selectors_ui <- renderUI({
    req(data_raw(), input$n_withhold)
    sel <- selections_r()
    cols <- setdiff(names(data_raw()), c(sel$sample, sel$batch, sel$class, sel$order))
    ids <- withheld_ids_r()
    if (!length(ids)) return(NULL)
    
    # snapshot previous selections **without** creating a dependency
    prev_all <- isolate(vapply(ids, function(id) input[[id]] %||% "", character(1)))
    
    lapply(seq_along(ids), function(i) {
      id   <- ids[i]
      prev <- prev_all[i]
      # prevent duplicates: exclude other picks, but keep this input's current value
      other_picks <- setdiff(prev_all, prev)
      choices_i <- c("Select a column..." = "", setdiff(cols, other_picks))
      
      selectInput(
        inputId  = id,
        label    = paste("Select column to withhold #", i),
        choices  = choices_i,
        selected = if (nzchar(prev) && prev %in% choices_i) prev else ""
      )
    })
  })
  
  #–– Cleaned data
  cleaned_r <- reactive({
    df  <- req(data_raw())
    sel <- selections_r()
    withheld <- withheld_r()
    
    # require 4 unique and  non-empty columns
    req(all(nzchar(c(sel$sample, sel$batch, sel$class, sel$order))))
    req(length(unique(c(sel$sample, sel$batch, sel$class, sel$order))) == 4)
    
    cleanData(df, sel$sample, sel$batch, sel$class, sel$order, withheld)
  }) %>% bindCache(input$file1$datapath, selections_r(), withheld_r())
  
  #–– Display Basic Information about data.
  output$basic_info <- renderUI({
    cleaned_data <- cleaned_r()
    req(cleaned_data)
    basicInfoUI(cleaned_data$df, cleaned_data$replacement_counts)
  })
  
  #-- Select a control class for fold change comparisons in corrected data.
  output$control_class_selector <- renderUI({
    req(cleaned_r())
    df <- cleaned_r()$df
    classes <- unique(df$class[df$class != "QC"])
    dropdown_choices <- c("Select a class..." = "", classes)
    
    tooltip(
      selectInput(
        "control_class",
        "Control Group",
        choices = dropdown_choices,
        selected = ""
      ),
      "Name of control samples in class column. This class's average will be used to compute fold changes in the corrected data file.",
      placement = "right"
    )
  })
  
  #–– Filter Data based on missing values.
  filtered_r <- reactive({
    cleaned_data <- req(cleaned_r())
    filter_data(cleaned_data$df, setdiff(
      names(cleaned_data$df),
      c("sample", "batch", "class", "order")
    ), input$Frule)
  }) #%>% bindCache(input$file1$datapath, selections_r(), input$Frule)
  output$filter_info <- renderUI({
    filtered_data <- filtered_r()
    req(filtered_data)
    filterInfoUI(filtered_data$mv_removed_cols)
  })
  
  #-- Move to next panel after inspecting the raw data
  observeEvent(input$next_correction, {
    req(cleaned_r(), filtered_r())
    
    # snapshot cleaned, filtered, and user inputs from Tab 1:
    rv$cleaned <- cleaned_r()
    rv$filtered <- filtered_r()
    sel <- selections_r()
    rv$params <- list(
      sample_col    = sel$sample,
      batch_col     = sel$batch,
      class_col     = sel$class,
      order_col     = sel$order,
      withheld_cols = sel$withheld,
      Frule         = input$Frule,
      control_class = input$control_class %||% NULL
    )
    
    updateTabsetPanel(session, "main_steps", "2. Correction Settings")
  })
  
  ########################## TAB 2
  
  #-- QC missing value warning and impute options.
  output$qc_missing_value_warning <- renderUI({
    req(rv$filtered)
    qcMissingValueWarning(rv$filtered$df)
  })
  output$qcImpute <- renderUI({
    req(rv$filtered)
    metab_cols <- setdiff(names(rv$filtered$df),
                          c('sample', 'batch', 'class', 'order'))
    qcImputeUI(rv$filtered$df, metab_cols)
  })
  
  #-- Sample missing value impute options.
  output$sampleImpute <- renderUI({
    req(rv$filtered)
    metab_cols <- setdiff(names(rv$filtered$df),
                          c('sample', 'batch', 'class', 'order'))
    sampleImputeUI(rv$filtered$df, metab_cols)
  })
  
  #-- Select correction method based on whats available for the data.
  output$correctionMethod <- renderUI({
    req(rv$filtered)
    correctionMethodUI(rv$filtered$df)
  })
  
  #-- Display unavailable options
  output$unavailable_options <- renderUI({
    req(rv$filtered)
    metab_cols <- setdiff(names(rv$filtered$df),
                          c('sample', 'batch', 'class', 'order'))
    unavailableOptionsUI(rv$filtered$df, metab_cols)
  })
  
  #-- reactive value for metab cols
  metab_cols_r <- reactive({
    req(rv$filtered)
    setdiff(names(rv$filtered$df), c("sample","batch","class","order"))
  })
  
  #-- reactives for NA value check
  has_qc_na_r <- reactive({
    req(rv$filtered)
    df <- rv$filtered$df
    mc <- metab_cols_r()
    any(is.na(dplyr::filter(df, .data$class == "QC")[, mc, drop = FALSE]))
  }) #%>% bindCache(input$file1$datapath, input$Frule)
  has_sam_na_r <- reactive({
    req(rv$filtered)
    df <- rv$filtered$df
    mc <- metab_cols_r()
    any(is.na(dplyr::filter(df, .data$class != "QC")[, mc, drop = FALSE]))
  }) #%>% bindCache(input$file1$datapath, input$Frule)
  
  #-- Impute missing values
  imputed_r <- reactive({
    req(rv$filtered)
    mc <- metab_cols_r()
    
    qcImpute  <- if (isTRUE(!has_qc_na_r()))  "nothing_to_impute" else input$qcImputeM
    samImpute <- if (isTRUE(!has_sam_na_r())) "nothing_to_impute" else input$samImputeM
    
    impute_missing(rv$filtered$df, mc, qcImpute, samImpute)
  }) #%>% bindCache(rv$filtered$df, input$qcImputeM, input$samImputeM) %>% identity()
  
  #-- Corrected data
  corrected_r <- eventReactive(input$correct, {
    req(imputed_r(), input$corMethod)
    
    # Isolate the current imputed result & inputs to "freeze" the snapshot at click time
    imputed <- isolate(imputed_r())
    cor_method <- isolate(input$corMethod)
    mc <- isolate(metab_cols_r())
    
    # take snapshot of imputed
    rv$imputed <- imputed
    rv$corrected <- correct_data(imputed$df, 
                              mc,
                              cor_method)
  })
  
  observeEvent(input$correct, { corrected_r()})
  
  #-- Filter corrected data
  filtered_corrected_r <- reactive({
    req(rv$filtered, rv$corrected)
    df_corrected <- rv$corrected$df
    
    if (isTRUE(input$remove_imputed)) {
      df <- remove_imputed_from_corrected(rv$filtered$df, df_corrected)
    } else {
      df <- df_corrected
    }
    
    if (input$post_cor_filter == FALSE) {
      rsd_filter(df,
                 input$rsd_filter,
                 c("sample", "batch", "class", "order"))
    } else {
      rsd_filter(df, Inf, c("sample", "batch", "class", "order"))
    }
  })
  
  #-- Option to withhold columns from TRN
  observe({
    req(rv$corrected, input$trn_withhold_checkbox)
    
    max_withhold <- max(ncol(rv$corrected$df) - 4, 0)
    
    output$trn_withhold_ui <- renderUI({
      if (input$transform == "TRN") {
        numericInput(
          inputId = "trn_withhold_n",
          label = "Number of columns to withold from TRN",
          value = 1,
          min = 1,
          max = max_withhold
        )
      }
    })
  })
  output$trn_withhold_selectors_ui <- renderUI({
    req(rv$corrected, input$trn_withhold_n)
    cols <- names(rv$corrected$df)
    cols <- setdiff(cols, c("sample", "batch", "class", "order"))
    dropdown_choices <- c("Select a column..." = "", cols)
    
    # Generate list of columns to withhold
    lapply(seq_len(input$trn_withhold_n), function(i) {
      selectInput(
        inputId = paste0("trn_withhold_col_", i),
        label = paste("Select column to withhold #", i),
        choices = dropdown_choices,
        selected = ""
      )
    })
  })
  
  #-- Scale/Normalize corrected data.
  transformed_r <- reactive({
    req(filtered_corrected_r())
    withheld_cols <- character(0)
    
    if (isTRUE(input$trn_withhold_checkbox) &&
        !is.null(input$trn_withhold_n)) {
      for (i in seq_len(input$trn_withhold_n)) {
        col <- input[[paste0("trn_withhold_col_", i)]]
        if (!is.null(col) &&
            col %in% names(filtered_corrected_r()$df)) {
          withheld_cols <- c(withheld_cols, col)
        }
      }
    }
    transform_data(filtered_corrected_r()$df,
                   input$transform,
                   withheld_cols,
                   input$ex_ISTD)
  })
  
  #-- Display corrected/transformed data and information.
  output$post_cor_filter_info <- renderUI({
    req(filtered_corrected_r())
    postCorFilterInfoUI(filtered_corrected_r())
  })
  output$cor_data <- renderTable({
    req(transformed_r())
    transformed_r()$df
  })
  
  #-- Download corrected data file only.
  output$download_corr_btn <- renderUI({
    req(transformed_r())
    
    div(
      style = "max-width: 300px; display: inline-block;",
      downloadButton(
        outputId = "download_corr_data",
        label    = "Download Excel File",
        class    = "btn btn-primary"
      )
    )
  })
  output$download_corr_data <- downloadHandler(
    filename = function() {
      paste0("corrected_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      fc <- isolate(filtered_corrected_r())
      tr <- isolate(transformed_r())
      
      wb <- corrected_file_download(
        input,
        rv$cleaned,
        rv$filtered,
        rv$imputed,
        rv$corrected,
        fc,
        tr
      )
      # Save to file
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  #-- Move to next tab after inspecting the corrected data
  observeEvent(input$next_visualization, {
    req(filtered_corrected_r(), transformed_r())
    # snapshot filtered corrected and transformed before making figures.
    rv$filtered_corrected <- isolate(filtered_corrected_r())
    rv$transformed <- isolate(transformed_r())
    
    updateTabsetPanel(session,
                      "main_steps",
                      "3. Evaluation Metrics and Visualization")
  })
  
  #-- display RSD comparison plot
  output$rsd_comparison_plot <- renderPlot(execOnResize = FALSE, res = 120,{
    print("RSD_COMPARISON_PLOT_ENTERED")
    req(rv$filtered, rv$filtered_corrected, input$rsd_cal)
    
    print("REQ CALL")
    
    df_before <- rv$filtered$df
    df_after <- rv$filtered_corrected$df
    rsd_mode  <- input$rsd_cal
    
    
    # Need at least 1 metabolite column
    validate(
      need(ncol(df_before) > 4L, "No metabolites left before correction."),
      need(ncol(df_after)  > 4L, "No metabolites left after correction."),
      need(sum(df_before$class == "QC", na.rm = TRUE) >= 2, "Not enough QC samples before correction (need >= 2)."),
      need(sum(df_after$class  == "QC",  na.rm = TRUE) >= 2, "Not enough QC samples after correction (need >= 2).")
    )
    
    print("POST-VALIDATE CALL 1")
    
    #grDevices::dev.hold(); on.exit(grDevices::dev.flush(), add = TRUE)
    
    #if ("package:gridExtra" %in% search()) detach("package:gridExtra", unload = TRUE, character.only = TRUE)
    
    # 3) Run all heavy work in isolate so any writes inside helpers cannot create deps.
    isolate({
      if (identical(rsd_mode, "met")) {
        p <- plot_rsd_comparison(df_before, df_after)
        print(paste("plot made:", class(p)))
      } else {
        p <- plot_rsd_comparison_class_met(df_before, df_after)
        print(paste("other plot made", class(p)))
      }
    })
    
    print(p)
  })
  
  #-- PCA plot
  output$pca_plot <- renderPlot({
    req(rv$imputed, rv$filtered_corrected)
    df <- rv$filtered_corrected$df
    mets <- setdiff(names(df), c("sample","batch","class","order"))
    validate(
      need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
      need(nrow(df) >= 3, "Need at least 3 samples for PCA.")
    )
    
    # Also ensure non-constant / non-NA columns
    X <- df[, mets, drop = FALSE]
    keep <- vapply(X, function(v) {
      v <- suppressWarnings(as.numeric(v))
      ok <- all(is.finite(v))
      nz <- (length(unique(v)) >= 2)
      ok && nz
    }, logical(1))
    validate(need(any(keep), "All metabolite columns are constant/invalid after filtering."))
    
    before <- rv$imputed
    after <- rv$filtered_corrected
    
    tryCatch({
      plot_pca(input, before, after, input$color_col)
    }, error = function(e) {
      showNotification(paste("PCA failed:", e$message),
                       type = "error", duration = 8)
      ggplot2::ggplot() + ggplot2::labs(title = "PCA failed — see notification")
    })
  }, res = 120)
  
  #-- Let user select which metabolite to display in scatter plot
  output$met_plot_selectors <- renderUI({
    req(rv$filtered, rv$transformed)
    raw_cols <- setdiff(names(rv$filtered$df),    c("sample","batch","class","order"))
    cor_cols <- setdiff(names(rv$transformed$df), c("sample","batch","class","order"))
    cols <- intersect(raw_cols, cor_cols)
    validate(need(length(cols) >= 1, "No overlapping metabolites between raw and corrected data."))
    selectInput("met_col", "Metabolite column", choices = cols, selected = cols[1])
  })
  
  output$metab_scatter <- renderPlot({
    req(input$met_col, rv$filtered, rv$transformed, rv$corrected)
    f_df <- rv$filtered$df
    t_df <- rv$transformed$df
    cor_method <- rv$corrected$str
    tryCatch({
      if (cor_method %in% c("Random Forest","Batchwise Random Forest")) {
        met_scatter_rf(f_df, t_df, i = input$met_col)
      } else if (cor_method %in% c("LOESS","Batchwise LOESS")) {
        met_scatter_loess(f_df, t_df, i = input$met_col)
      } else {
        ggplot2::ggplot() + ggplot2::labs(title = "No correction method selected.")
      }
    }, error = function(e) {
      showNotification(paste("Scatter failed:", e$message),
                       type = "error", duration = 8)
      ggplot2::ggplot() + ggplot2::labs(title = "Scatter failed — see notification")
    })
  })
  
  #-- Download all figures as zip folder.
  output$download_fig_zip_btn <- renderUI({
    req(rv$transformed)
    
    div(
      style = "max-width: 300px; display: inline-block;",
      downloadButton("download_fig_zip", "Download All Figures", class = "btn-primary")
    )
  })
  # -- progress bar
  progress_reactive <- reactiveVal(0)
  #-- progress for downloading all images
  output$progress_ui <- renderUI({
    req(progress_reactive() > 0, progress_reactive() <= 1)
    
    div(
      style = "margin-top: 10px;",
      tags$label("Progress:"),
      tags$progress(
        value = progress_reactive(),
        max = 1,
        style = "width: 100%; height: 20px;"
      ),
      tags$span(sprintf("%.0f%%", progress_reactive() * 100))
    )
  })
  
  output$download_fig_zip <- downloadHandler(
    filename = function() {
      paste0("figures_", Sys.Date(), ".zip")
    },
    content = function(file) {
      figs <- figure_folder_download(input, rv$imputed, rv$filtered, rv$filtered_corrected)
      
      zipfile <- tempfile(fileext = ".zip")
      old_wd <- setwd(figs$tmp_dir)
      on.exit({
        unlink(figs$fig_dir, recursive = TRUE)
        unlink(zipfile)
        setwd(old_wd)
      }, add = TRUE)
      zip(zipfile = zipfile,
          files = "figures",
          extras = "-r9Xq")
      file.copy(zipfile, file)
      
      # Remove progress bar
      progress_reactive(0)
    }
  )
  
  #-- Move to next tab after inspecting the corrected data figures
  observeEvent(input$next_export, {
    updateTabsetPanel(session, "main_steps", "4. Export Corrected Data and Plots")
  })
  
  output$download_all_ui <- renderUI({
    req(rv$transformed)
    downloadButton(
      outputId = "download_all_zip",
      label = "Download All",
      class = "btn-primary btn-lg"
    )
  })
  
  #-- Allow user to download corrected data, figures, and correction report.
  output$download_all_zip <- downloadHandler(
    filename = function() {
      paste0("corrected_data_plots_report_", Sys.Date(), ".zip")
    },
    content = function(file) {
      base_dir <- tempfile()
      dir.create(base_dir)
      
      # 1. create and save corrected data file
      cor_data_filename <- paste0("corrected_data_", Sys.Date(), ".xlsx")
      cor_data_path <- file.path(base_dir, cor_data_filename)
      wb <- corrected_file_download(
        input,
        rv$cleaned,
        rv$filtered,
        rv$imputed,
        rv$corrected,
        rv$filtered_corrected,
        rv$transformed
      )
      saveWorkbook(wb, cor_data_path, overwrite = TRUE)
      
      # 2. create and save figure folder
      fig_info <- figure_folder_download(input, rv$imputed, rv$filtered, rv$filtered_corrected)
      fig_dir <- fig_info$fig_dir
      figures_path <- file.path(base_dir, "figures")
      dir.create(figures_path)
      file.copy(from = list.files(fig_dir, full.names = TRUE),
                to = figures_path,
                recursive = TRUE)
      
      # make zip file
      zipfile <- tempfile(fileext = ".zip")
      old_wd <- setwd(base_dir)
      on.exit({
        unlink(base_dir, recursive = TRUE)
        unlink(zipfile)
        setwd(old_wd)
      }, add = TRUE)
      zip(zipfile = zipfile,
          files = list.files(base_dir),
          extras = "-r9Xq")
      
      file.copy(zipfile, file)
    }
  )
  
}


shinyApp(ui = ui, server = server)
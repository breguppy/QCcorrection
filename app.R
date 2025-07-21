library(shiny)
library(shinyBS)
library(bslib)
library(dplyr)
library(shinycssloaders)
library(openxlsx)
source("R/helpers.R")
source("R/processing_helpers.R")
source("R/met_scatter_rf.R")
source("R/met_scatter_loess.R")
source("R/download_helpers.R")


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
               tags$h6(
                 "When sorted by injection order, the data must begin and end with a QC samples."
               ),
               fileInput(
                 inputId = "file1",
                 label = "Choose CSV File",
                 accept = ".csv",
                 buttonLabel = "Browse...",
                 placeholder = "No file selected"
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
      card(
        style = "background-color: #eeeeee;",
        tags$h4("1.3 Identify Control Group"),
        fluidRow(
          column(4, tooltip(
            checkboxInput(
              inputId = "no_control",
              label = "No control group.",
              value = FALSE
            ),
            "Check the box if The data does not have a control group.",
            placement = "right"
          )),
          column(4, conditionalPanel(
            "input.no_control == false",
            uiOutput("control_class_selector")
          ))
        )
      ),
      #-- Filter data based on missing values
      card(layout_sidebar(
        sidebar = sidebar(
          #--- Filter metabolites
          tags$h4("1.4 Filter Raw Data"),
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
        uiOutput("filter_info")
      ), ),
      #-- move on to step 2
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
          column(
            3,
            tags$h5("Impute Missing QC Values"),
            uiOutput("qcImpute")
          ),
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
                "Total Ratio Normalization (TRN)" = "TRN",
                "None" = "none"
              ),
              selected = "none"
            ),
            conditionalPanel(
              "input.transform == 'TRN'",
              tooltip(
                checkboxInput(
                  inputId = "ex_ISTD",
                  label = "Exclude Internal Standards from TRN.",
                  value = TRUE
                ),
                "Check this box if you do not want internal standards to be included in the TRN calculation.",
                placement = "right"
              )
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
        ),
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
        plotOutput("rsd_comparison_plot", height = "600px", width = "1100px")
      )),
      card(layout_sidebar(
        sidebar = sidebar(
          #--- plot metabolite
          tags$h4("3.2 Metabolite Scatter plots"),
          uiOutput("met_plot_selectors"),
          width = 400,
        ),
        plotOutput("metab_scatter", height = "600px", width = "600px"),
      )),
      card(
        layout_sidebar(
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
          uiOutput("download_fig_zip_btn", container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          uiOutput("progress_ui"),
        )
      ),
      card(
        actionButton(
          inputId = "next_export",
          label = "Next: Export Corrected Data and Plots",
          class = "btn-primary btn-lg"
        ),
      )
    ),
    #--- Step 5: Export Data
    nav_panel(
      title = "4. Export Corrected Data and Plots", 
      
      card(
        card_title("Download Data and Plots"),
        tags$span("TODO: Describe what will be downloaded and the format."),
        downloadButton(outputId = "download_all",
                       label = "Download All", 
                       class = "btn-primary btn-lg"),
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  #-- Import Raw Data
  data_raw <- reactive({
    req(input$file1)
    #colnames_original <- names(read.csv(input$file1$datapath, nrows = 1, check.names = FALSE))
    df <- read.csv(input$file1$datapath, header = TRUE, check.names = FALSE)
  })
  
  #–– View Raw Data
  output$contents <- renderTable(data_raw())
  
  #–– Column selection & warnings
  output$column_selectors <- renderUI({
    req(data_raw())
    cols <- names(data_raw())
    
    dropdown_choices <- c("Select a column..." = "", cols)
    
    tagList(
      tooltip(
        selectInput(
          "sample_col",
          "sample column",
          choices = dropdown_choices,
          selected = ""
        ),
        "Column that contains unique sample names.",
        placement = "right"
      ),
      tooltip(
        selectInput(
          "batch_col",
          "batch column",
          choices = dropdown_choices,
          selected = ""
        ),
        "Column that contains batch information.",
        placement = "right"
      ),
      tooltip(
        selectInput(
          "class_col",
          "class column",
          choices = dropdown_choices,
          selected = ""
        ),
        "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.",
        placement = "right"
      ),
      tooltip(
        selectInput(
          "order_col",
          "order column",
          choices = dropdown_choices,
          selected = ""
        ),
        "Column that indicates injection order.",
        placement = "right"
      )
    )
  })
  output$column_warning <- renderUI({
    req(data_raw())
    selected <- c(input$sample_col,
                  input$batch_col,
                  input$class_col,
                  input$order_col)
    columnWarningUI(data_raw(), selected)
  })
  
  #-- Option to withhold more columns form correction.
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
    req(data_raw(), input$withhold_cols, input$n_withhold)
    cols <- names(data_raw())
    cols <- setdiff(cols,
                    c(
                      input$sample_col,
                      input$batch_col,
                      input$class_col,
                      input$order_col
                    ))
    dropdown_choices <- c("Select a column..." = "", cols)
    
    # Generate list of columns to withhold from correction.
    lapply(seq_len(input$n_withhold), function(i) {
      selectInput(
        inputId = paste0("withhold_col_", i),
        label = paste("Select column to withhold #", i),
        choices = dropdown_choices,
        selected = ""
      )
    })
  })
  
  #–– Cleaned data
  cleaned <- reactive({
    sel <- c(input$sample_col,
             input$batch_col,
             input$class_col,
             input$order_col)
    
    req(!any(sel == ""), length(unique(sel)) == 4)
    
    df <- data_raw()
    if (isTRUE(input$withhold_cols) && !is.null(input$n_withhold)) {
      withheld_cols <- character(0)
      for (i in seq_len(input$n_withhold)) {
        col <- input[[paste0("withhold_col_", i)]]
        if (!is.null(col) && col %in% names(df)) {
          withheld_cols <- c(withheld_cols, col)
        }
      }
      df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]
    }
    cleanData(df, sel[1], sel[2], sel[3], sel[4])
  })
  
  #–– Display Basic Information about data.
  output$basic_info <- renderUI({
    cleaned_data <- cleaned()
    req(cleaned_data)
    basicInfoUI(cleaned_data$df, cleaned_data$replacement_counts)
  })
  
 #-- Select a control class for fold change comparisons in corrected data.
  output$control_class_selector <- renderUI({
    req(cleaned())
    df <- cleaned()$df
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
  filtered <- reactive({
    cleaned_data <- cleaned()
    req(cleaned_data)
    filter_data(cleaned_data$df, setdiff(
      names(cleaned_data$df),
      c("sample", "batch", "class", "order")
    ), input$Frule)
  })
  output$filter_info <- renderUI({
    filtered_data <- filtered()
    req(filtered_data)
    filterInfoUI(filtered_data$removed_cols)
  })
  
  #-- Move to next panel after inspecting the raw data
  observeEvent(input$next_correction, {
    updateTabsetPanel(session, "main_steps", "2. Correction Settings")
  })
  
  #-- QC missing value warning and impute options.
  output$qc_missing_value_warning <- renderUI({
    req(filtered())
    qcMissingValueWarning(filtered()$df)
  })
  output$qcImpute <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    qcImputeUI(filtered()$df, metab_cols)
  })
  
  #-- Sample missing value impute options.
  output$sampleImpute <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    sampleImputeUI(filtered()$df, metab_cols)
  })
  
  #-- Select correction method based on whats available for the data.
  output$correctionMethod <- renderUI({
    req(filtered())
    correctionMethodUI(filtered()$df)
  })
  
  #-- Display unavailable options
  output$unavailable_options <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    unavailableOptionsUI(filtered()$df, metab_cols)
  })
  
  #-- Impute missing values
  imputed <- reactive({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    qc_df <- filtered()$df %>% filter(filtered()$df$class == "QC")
    has_qc_na <- any(is.na(qc_df[, metab_cols]))
    sam_df <- filtered()$df %>% filter(filtered()$df$class != "QC")
    has_sam_na <- any(is.na(sam_df[, metab_cols]))
    
    ifelse(!has_qc_na, qcImpute <- "nothing_to_impute", qcImpute <- input$qcImputeM)
    ifelse(!has_sam_na, samImpute <- "nothing_to_impute", samImpute <- input$samImputeM)
    
    impute_missing(filtered()$df,
                   setdiff(
                     names(filtered()$df),
                     c("sample", "batch", "class", "order")
                   ),
                   qcImpute, samImpute)
  })
  
  #-- Corrected data 
  corrected <- eventReactive(input$correct, {
    req(imputed())
    
    correct_data(imputed()$df,
                 setdiff(
                   names(imputed()$df),
                   c("sample", "batch", "class", "order")
                 ),
                 input$corMethod)
  })
  
  #-- Filter corrected data
  filtered_corrected <- reactive({
    req(filtered(), corrected())
    df_corrected <- corrected()$df
    
    if (isTRUE(input$remove_imputed)) {
      df <- remove_imputed_from_corrected(filtered()$df, df_corrected)
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
    req(corrected(), input$trn_withhold_checkbox)
    
    max_withhold <- max(ncol(corrected()$df) - 4, 0)
    
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
    req(corrected(), input$trn_withhold_n)
    cols <- names(corrected()$df)
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
  transformed <- reactive({
    req(filtered_corrected())
    withheld_cols <- character(0)
    
    if (isTRUE(input$trn_withhold_checkbox) && !is.null(input$trn_withhold_n)) {
      for (i in seq_len(input$trn_withhold_n)) {
        col <- input[[paste0("trn_withhold_col_", i)]]
        if (!is.null(col) && col %in% names(filtered_corrected()$df)) {
          withheld_cols <- c(withheld_cols, col)
        }
      }
    }
    transform_data(filtered_corrected()$df, input$transform, withheld_cols, input$ex_ISTD)
  })
  
  #-- Display corrected/transformed data and information.
  output$post_cor_filter_info <- renderUI({
    req(filtered_corrected())
    postCorFilterInfoUI(filtered_corrected())
  })
  output$cor_data <- renderTable({
    req(transformed())
    transformed()$df
  })
  
  #-- Download corrected data file only.
  output$download_corr_btn <- renderUI({
    req(transformed())
    
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
      wb <- corrected_file_download(input, 
                              filtered(), 
                              imputed(), 
                              corrected(), 
                              filtered_corrected(), 
                              transformed())
      # Save to file
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  #-- Move to next tab after inspecting the corrected data
  observeEvent(input$next_visualization, {
    updateTabsetPanel(session,
                      "main_steps",
                      "3. Evaluation Metrics and Visualization")
  })
  
  #-- display RSD comparison plot
  output$rsd_comparison_plot <- renderPlot({
    req(filtered(), filtered_corrected(), input$rsd_cal)
    if (input$rsd_cal == "met") {
      plot_rsd_comparison(filtered()$df,
                          filtered_corrected()$df)
    } else if (input$rsd_cal == "class_met") {
      plot_rsd_comparison_class_met(filtered()$df,
                                    filtered_corrected()$df)
    }
  })
  
  #-- Let user select which metabolite to display in scatter plot
  output$met_plot_selectors <- renderUI({
    req(filtered(), transformed())
    raw_cols <- setdiff(names(filtered()$df),
                        c("sample", "batch", "class", "order"))
    cor_cols <- setdiff(names(transformed()$df),
                        c("sample", "batch", "class", "order"))
    cols <- intersect(raw_cols, cor_cols)
    tagList(selectInput(
      "met_col",
      "Metabolite column",
      choices = cols,
      selected = cols[1]
    ))
  })
  output$metab_scatter <- renderPlot({
    req(input$met_col, filtered(), transformed(), input$corMethod)
    if (input$corMethod %in% c("RF", "BW_RF")) {
      met_scatter_rf(
        data_raw = filtered()$df,
        data_cor = transformed()$df,
        i = input$met_col
      )
    } else if (input$corMethod %in% c("LOESS", "BW_LOESS")) {
      met_scatter_loess(
        data_raw = filtered()$df,
        data_cor = transformed()$df,
        i = input$met_col
      )
    }
  })
  
  #-- Download all figures as zip folder.
  output$download_fig_zip_btn <- renderUI({
    req(transformed())
    
    div(
      style = "max-width: 300px; display: inline-block;",
      downloadButton("download_fig_zip", 
                     "Download All Figures", 
                     class = "btn-primary")
    )
  })
  # -- progress bar
  progress_reactive <- reactiveVal(0)
  #-- progress for downloading all images
  output$progress_ui <- renderUI({
    req(progress_reactive() > 0, progress_reactive() <= 1)
    
    tags$div(
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
      figs <- figure_folder_download(input, filtered(), filtered_corrected())
      
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
    updateTabsetPanel(session,
                      "main_steps",
                      "4. Export Corrected Data and Plots")
  })
  
  #-- Allow user to download corrected data, figures, and correction report.
  
}


shinyApp(ui = ui, server = server)
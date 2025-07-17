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
  data_raw <- reactive({
    req(input$file1)
    #colnames_original <- names(read.csv(input$file1$datapath, nrows = 1, check.names = FALSE))
    df <- read.csv(input$file1$datapath, header = TRUE, check.names = FALSE)
  })
  
  #––  preview
  output$contents <- renderTable(data_raw())
  
  #–– column selection & warning
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
    
    # Generate list of selectInputs
    lapply(seq_len(input$n_withhold), function(i) {
      selectInput(
        inputId = paste0("withhold_col_", i),
        label = paste("Select column to withhold #", i),
        choices = dropdown_choices,
        selected = ""
      )
    })
  })
  
  #–– cleaned data
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
  
  #–– basic info
  output$basic_info <- renderUI({
    cleaned_data <- cleaned()
    req(cleaned_data)
    basicInfoUI(cleaned_data$df, cleaned_data$replacement_counts)
  })
  
 
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
      "Name of control samples in class column.",
      placement = "right"
    )
  })
  
  #–– filter step
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
  output$qc_missing_value_warning <- renderUI({
    req(filtered())
    qcMissingValueWarning(filtered()$df)
  })
  
  output$qcImpute <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    qcImputeUI(filtered()$df, metab_cols)
  })
  
  output$sampleImpute <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    sampleImputeUI(filtered()$df, metab_cols)
  })
  
  output$correctionMethod <- renderUI({
    req(filtered())
    correctionMethodUI(filtered()$df)
  })
  
  output$unavailable_options <- renderUI({
    req(filtered())
    metab_cols <- setdiff(names(filtered()$df), c('sample', 'batch', 'class', 'order'))
    unavailableOptionsUI(filtered()$df, metab_cols)
  })
  
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
  
  corrected <- eventReactive(input$correct, {
    req(imputed())
    
    correct_data(imputed()$df,
                 setdiff(
                   names(imputed()$df),
                   c("sample", "batch", "class", "order")
                 ),
                 input$corMethod)
  })
  
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
    
    # Generate list of selectInputs
    lapply(seq_len(input$trn_withhold_n), function(i) {
      selectInput(
        inputId = paste0("trn_withhold_col_", i),
        label = paste("Select column to withhold #", i),
        choices = dropdown_choices,
        selected = ""
      )
    })
  })
  
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
  
  output$post_cor_filter_info <- renderUI({
    req(filtered_corrected())
    postCorFilterInfoUI(filtered_corrected())
  })
  
  output$cor_data <- renderTable({
    req(transformed())
    transformed()$df
  })
  
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
      # Create a new workbook
      wb <- createWorkbook()
      # make column names bold
      bold_style <- createStyle(textDecoration = "Bold")
      
      # Add 0. Raw Data tab
      addWorksheet(wb, "0. Raw Data")
      writeData(wb, sheet = "0. Raw Data", x = filtered()$df, startRow = 3, headerStyle = bold_style)
      
      # Add 1. Correction Settings
      correction_settings_df <- data.frame(
        Settings = c(
          "Sample Column Name",
          "Batch Column Name",
          "Class Column Name",
          "Order Column Name",
          "Missing Value Threshold",
          "QC Missing Value Imputation Method",
          "Sample Missing Value Imputation Method",
          "Correction Method",
          "Remove Imputed Values After Correction?",
          "QC RSD% Threshold",
          "Scaling/Tranformation Method"
        ),
        Values = c(
          input$sample_col,
          input$batch_col,
          input$class_col,
          input$order_col,
          paste0(filtered()$Frule, "%"),
          imputed()$qc_str,
          imputed()$sam_str,
          corrected()$str,
          input$remove_imputed,
          paste0(filtered_corrected()$rsd_cutoff, "%"),
          transformed()$str
        ),
        stringsAsFactors = FALSE
      )
      addWorksheet(wb, "1. Correction Settings")
      # TODO: add description at top of tab.
      writeData(wb, sheet = "1. Correction Settings", x = correction_settings_df, startRow = 3, startCol = 1, headerStyle = bold_style)
      
      # Append Missing Value Filtered Metabolites
      if (length(filtered()$removed_cols) > 0) {
        mv_df <- data.frame(
          Missing_Value_Filtered_Metabolites = filtered()$removed_cols,
          stringsAsFactors = FALSE
        )
        writeData(wb, "1. Correction Settings", x = mv_df, startRow = 3, startCol = 4, headerStyle = bold_style)
      }
      
      # Append QC RSD Filtered Metabolites
      if (length(filtered_corrected()$removed_metabolites) > 0) {
        rsd_df <- data.frame(
          QC_RSD_Filtered_Metabolites = filtered_corrected()$removed_metabolites,
          stringsAsFactors = FALSE
        )
        writeData(wb, "1. Correction Settings", x = rsd_df, startRow = 3, startCol = 6, headerStyle = bold_style)
      }
      
      if (transformed()$str == "TRN" && length(transformed()$withheld_cols) > 0) {
        ex_trn <- data.frame(
          Exculded_In_TRN = transformed()$withheld_cols,
          stringsAsFactors = FALSE
        )
        writeData(wb, "1. Correction Settings", x = ex_trn, startRow = 3, startCol = 8, headerStyle = bold_style)
      }
      
      # Add 2. Drift Normalized tab
      corrected_df <- filtered_corrected()$df
      samples <- corrected_df[corrected_df$class != "QC", ]
      addWorksheet(wb, "2. Drift Normalized")
      # TODO: add description at top of tab.
      writeData(wb, sheet = "2. Drift Normalized", x = samples, startRow = 3, headerStyle = bold_style)
      
      # Add Scaled Data tab (name depends on method)
      transformed_df <- transformed()$df
      keep_cols <- setdiff(names(transformed_df), transformed()$withheld_cols)
      transformed_df <- transformed_df[transformed_df$class != "QC" , keep_cols]
      addWorksheet(wb, "3. Scaled or Normalized")
      # TODO: add description at top of tab.
      writeData(wb, sheet = "3. Scaled or Normalized", x = transformed_df, startRow = 3, headerStyle = bold_style)
      
      # Add 4. Grouped Data Organized
      grouped_data <- group_stats(transformed_df)
      addWorksheet(wb, "4. Grouped Data Organized")
      current_row <- 3
      for (group_name in names(grouped_data$group_dfs)) {
        group <- grouped_data$group_dfs[[group_name]]
        group_size <- nrow(group)
        
        writeData(wb, sheet = "4. Grouped Data Organized", x = group, startRow = current_row, headerStyle = bold_style)
        current_row <- current_row + group_size + 1
        group_stats <- grouped_data$group_stats_dfs[[group_name]]
        writeData(wb, sheet = "4. Grouped Data Organized", x = group_stats, startRow = current_row, startCol = 2, headerStyle = bold_style)
        current_row <- current_row + 6
      }
      
      # Add. 5. Grouped Data Fold Change
      control_stats <- grouped_data$group_stats_dfs[[input$control_class]]
      fold_change <- fold_changes(transformed_df, control_stats[1, ])
      group_fc_data <- group_stats(fold_change)
      addWorksheet(wb, "5. Group Data Fold Change")
      current_row <- 3
      for (group_name in names(group_fc_data$group_dfs)) {
        group <- group_fc_data$group_dfs[[group_name]]
        group_size <- nrow(group)
        
        writeData(wb, sheet = "5. Group Data Fold Change", x = group, startRow = current_row, headerStyle = bold_style)
        current_row <- current_row + group_size + 1
        group_stats <- group_fc_data$group_stats_dfs[[group_name]]
        writeData(wb, sheet = "5. Group Data Fold Change", x = group_stats, startRow = current_row, startCol = 2, headerStyle = bold_style)
        current_row <- current_row + 6
      }
      
      # Add. Appendix1. Metaboanalyst Ready
      names(fold_change)[names(fold_change) == "sample"] <- "Sample Name"
      names(fold_change)[names(fold_change) == "class"] <- "Group"
      fold_change$batch <- NULL
      fold_change$order <- NULL
      addWorksheet(wb, "Appendix1. Metaboanalyst Ready")
      writeData(wb, sheet = "Appendix1. Metaboanalyst Ready", x = fold_change)
       
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
      # create temp folder for figures
      tmp_dir <- tempdir()
      fig_dir <- file.path(tmp_dir, "figures")
      if (dir.exists(fig_dir))
        unlink(fig_dir, recursive = TRUE)
      dir.create(fig_dir)
      # create RSD figure folder
      rsd_fig_dir <- file.path(fig_dir, "RSD figures")
      if (dir.exists(rsd_fig_dir))
        unlink(rsd_fig_dir, recursive = TRUE)
      dir.create(rsd_fig_dir)
      
      # create metabolite figure folder
      met_fig_dir <- file.path(fig_dir, "metabolite figures")
      if (dir.exists(met_fig_dir))
        unlink(met_fig_dir, recursive = TRUE)
      dir.create(met_fig_dir)
      
      # create RSD plots
      if (input$rsd_cal == "met") {
        rsd_fig <- plot_rsd_comparison(filtered()$df,
                                       filtered_corrected()$df)
      } else if (input$rsd_cal == "class_met") {
        rsd_fig <- plot_rsd_comparison_class_met(filtered()$df,
                                                 filtered_corrected()$df)
      }
      rsd_path <- file.path(rsd_fig_dir,
                            paste0("rsd_comparison_", input$rsd_cal, ".", input$fig_format))
      if (input$fig_format == "png") {
        ggsave(
          rsd_path,
          plot = rsd_fig,
          width = 16,
          height = 8,
          dpi = 300
        )
      } else if (input$fig_format == "pdf") {
        ggsave(rsd_path,
               plot = rsd_fig,
               width = 16,
               height = 8)
      }
      
      # create metabolite scatter plots
      raw_cols <- setdiff(names(filtered()$df),
                          c("sample", "batch", "class", "order"))
      cor_cols <- setdiff(names(filtered_corrected()$df),
                          c("sample", "batch", "class", "order"))
      cols <- intersect(raw_cols, cor_cols)
      n <- length(cols)
      withProgress(message = "Creating figures...", value = 0, {
        for (i in seq_along(cols)) {
          metab <- cols[i]
          if (input$corMethod %in% c("RF", "BW_RF")) {
            fig <- met_scatter_rf(
              data_raw = filtered()$df,
              data_cor = filtered_corrected()$df,
              i = metab
            )
          } else if (input$corMethod %in% c("LOESS", "BW_LOESS")) {
            fig <- met_scatter_loess(
              data_raw = filtered()$df,
              data_cor = filtered_corrected()$df,
              i = metab
            )
          }
          metab <- sanitize_figname(metab)
          path <- file.path(met_fig_dir, paste0(metab, ".", input$fig_format))
          if (input$fig_format == "png") {
            ggsave(
              path,
              plot = fig,
              width = 8,
              height = 8,
              dpi = 300
            )
          } else if (input$fig_format == "pdf") {
            ggsave(path,
                   plot = fig,
                   width = 8,
                   height = 8)
          }
          incProgress(1 / n, detail = paste("Saved:", metab))
        }
      })
      
      # Zip the folder
      zipfile <- tempfile(fileext = ".zip")
      old_wd <- setwd(tmp_dir)
      on.exit({
        unlink(fig_dir, recursive = TRUE)
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
  observeEvent(input$next_export, {
    updateTabsetPanel(session,
                      "main_steps",
                      "4. Export Corrected Data and Plots")
  })
}


shinyApp(ui = ui, server = server)
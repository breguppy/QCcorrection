library(shiny)
library(shinyBS)
library(bslib)
library(dplyr)
library(shinycssloaders)
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
            "Metabolites with more then the acceptable % of missing values will be removed from the data.",
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
            radioButtons(
              inputId = "qcImputeM",
              label = "QC Imputation Method",
              choices = list(
                "metabolite median" = "median",
                "metabolite mean" = "mean",
                "class-metabolite median" = "class_median",
                "class-metabolite mean" = "class_mean",
                "minimum value" = "min",
                "half minimum value" = "minHalf",
                "KNN" = "KNN",
                "zero" = "zero"
                ),
            selected = "median",
            inline = FALSE
            )
          ),
          column(
            3,
            tags$h5("Impute Missing Sample Values"),
            radioButtons(
              inputId = "samimputeM",
              label = "Sample Imputation Method",
              choices = list(
                "metabolite median" = "median",
                "metabolite mean" = "mean",
                "QC-metabolite median" = "class_median",
                "QC-metabolite mean" = "class_mean",
                "minimum value" = "min",
                "half minimum value" = "minHalf",
                "KNN" = "KNN",
                "zero" = "zero"
              ),
              selected = "median",
              inline = FALSE
            )
          ),
          column(
            3,
            #--- Choose Correction method
            tags$h5("Choose Correction Method"),
            radioButtons(
              inputId = "corMethod",
              label = "Method",
              choices = list(
                "Random Forest" = "RF",
                "Local Polynomial Fit (LOESS)" = "LOESS",
                "Batchwise Random Forest" = "BW_RF",
                "Batchwise Local polynomial fit (LOESS)" = "BW_LOESS"
              ),
              selected = "RF"
            ),
          ),
          actionButton(
            inputId = "correct",
            label = "Correct Data with Selected Settings",
            class = "btn-primary btn-lg"
          ),
          #uiOutput("correction_info"),
        ),
      ),
      #-- User selects post-correction filtering and transformation/normalization
      card(
        layout_sidebar(
          sidebar = sidebar(
            # After correction filtering
            tags$h4("2.2 Post-Correction Filtering"),
            checkboxInput(
              inputId = "remove_imputed",
              label = "Remove imputed values after correction?",
              value = FALSE
            ),
            checkboxInput(
              inputId = "post_cor_filter",
              label = "Don't filter metabolites based on QC RSD%",
              value = FALSE
            ),
            conditionalPanel(
              "input.post_cor_filter == false",
              sliderInput(
                inputId = "rsd_filter",
                label = "Metabolite RSD% threshold for QC samples",
                min = 0,
                max = 100,
                value = 50
              )
            ),
            tags$hr(),
            
            # After correction scaling / normalization
            tags$h4("2.3 Post-Correction Transformation or Normalization"),
            tags$h6(style = "color: darkorange; font-weight: bold;", "(Coming Soon!)"),
            radioButtons(
              inputId = "transform",
              label = "Method",
              choices = list(
                "Log 2 transformation" = "log2",
                "Total Ratio Normalization" = "TRN",
                "None" = "none"
              ),
              selected = "none"
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
    
    #--- Step 4: plots
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
          downloadButton("download_fig_zip", "Download All Figures", class = "btn-primary btn-lg"),
          uiOutput("progress_ui"),
        )
      )
    ),
    #--- Step 5: Export Data
    nav_panel(title = "4. Export Corrected Data and Plots", card(card_title(
      "TODO: Export button"
    ), ))
  )
)

# Define server
server <- function(input, output, session) {
  data_raw <- reactive({
    req(input$file1)
    df <- read.csv(input$file1$datapath, header = TRUE)
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
    filtered_data <- filtered()
    req(filtered_data)
    qcMissingValueWarning(filtered_data$df_filtered)
  })
  imputed <- reactive({
    filtered_result <- filtered()
    req(filtered_result)
    
    impute_missing(filtered_result$df_filtered,
                   setdiff(
                     names(filtered_result$df_filtered),
                     c("sample", "batch", "class", "order")
                   ),
                   input$imputeM)
  })
  
  output$correction_info <- renderUI({
    imputed_result <- imputed()
    req(imputed_result)
    correctionInfoUI(imputed_result, input$imputeM, input$corMethod)
  })
  
  corrected <- eventReactive(input$correct, {
    imputed_result <- imputed()
    req(imputed_result)
    
    correct_data(imputed_result$df_imputed,
                 setdiff(
                   names(imputed_result$df_imputed),
                   c("sample", "batch", "class", "order")
                 ),
                 input$corMethod)
  })
  
  filtered_corrected <- reactive({
    df_corrected <- corrected()
    filtered_result <- filtered()
    req(df_corrected)
    if (isTRUE(input$remove_imputed)) {
      df <- remove_imputed_from_corrected(filtered_result$df_filtered, df_corrected)
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
    
    max_withhold <- max(ncol(corrected()) - 4, 0)
    
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
    cols <- names(corrected())
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
    fil_cor_result <- filtered_corrected()
    req(fil_cor_result)
    withheld_cols <- character(0)
    if (isTRUE(input$trn_withhold_checkbox) && !is.null(input$trn_withhold_n)) {
      for (i in seq_len(input$trn_withhold_n)) {
        col <- input[[paste0("trn_withhold_col_", i)]]
        if (!is.null(col) && col %in% names(df)) {
          withheld_cols <- c(withheld_cols, col)
        }
      }
    }
    transform_data(fil_cor_result$filtered_df, input$transform, withheld_cols)
  })
  
  output$post_cor_filter_info <- renderUI({
    fil_cor_result <- filtered_corrected()
    req(fil_cor_result)
    postCorFilterInfoUI(fil_cor_result)
  })
  output$cor_data <- renderTable({
    req(transformed())
    transformed()
  })
  
  output$download_corr_btn <- renderUI({
    req(transformed())
    
    downloadButton(outputId = "download_corr_data",
                   label    = "Download CSV",
                   class    = "btn btn-primary")
  })
  
  output$download_corr_data <- downloadHandler(
    filename = function() {
      paste0("corrected_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_corrected()$filtered_df, file, row.names = FALSE)
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
      plot_rsd_comparison(filtered()$df_filtered,
                          filtered_corrected()$filtered_df)
    } else if (input$rsd_cal == "class_met") {
      plot_rsd_comparison_class_met(filtered()$df_filtered,
                                    filtered_corrected()$filtered_df)
    }
  })
  
  #-- Let user select which metabolite to display in scatter plot
  output$met_plot_selectors <- renderUI({
    req(filtered(), transformed())
    raw_cols <- setdiff(names(filtered()$df_filtered),
                        c("sample", "batch", "class", "order"))
    cor_cols <- setdiff(names(transformed()),
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
        data_raw = filtered()$df_filtered,
        data_cor = transformed(),
        i = input$met_col
      )
    } else if (input$corMethod %in% c("LOESS", "BW_LOESS")) {
      met_scatter_loess(
        data_raw = filtered()$df_filtered,
        data_cor = transformed(),
        i = input$met_col
      )
    }
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
        rsd_fig <- plot_rsd_comparison(filtered()$df_filtered,
                                       filtered_corrected()$filtered_df)
      } else if (input$rsd_cal == "class_met") {
        rsd_fig <- plot_rsd_comparison_class_met(filtered()$df_filtered,
                                                 filtered_corrected()$filtered_df)
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
      raw_cols <- setdiff(names(filtered()$df_filtered),
                          c("sample", "batch", "class", "order"))
      cor_cols <- setdiff(names(filtered_corrected()$filtered_df),
                          c("sample", "batch", "class", "order"))
      cols <- intersect(raw_cols, cor_cols)
      n <- length(cols)
      withProgress(message = "Creating figures...", value = 0, {
        for (i in seq_along(cols)) {
          metab <- cols[i]
          if (input$corMethod == "RF") {
            fig <- met_scatter_rf(
              data_raw = filtered()$df_filtered,
              data_cor = filtered_corrected()$filtered_df,
              i = metab
            )
          } else if (input$corMethod == "LOESS") {
            fig <- met_scatter_loess(
              data_raw = filtered()$df_filtered,
              data_cor = filtered_corrected()$filtered_df,
              i = metab
            )
          }
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
}


shinyApp(ui = ui, server = server)
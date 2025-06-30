# Load packages

library(shiny)
library(bslib)
library(dplyr)
library(shinycssloaders)
source("helpers.R")
source("processing_helpers.R")
source("plotting_helpers.R")


# Define UI

ui <- fluidPage(
  theme = bs_theme(preset = "cosmo"),
  
  accordion(
    id = "main_steps",
    
    #--- Step 1: Raw Data
    accordion_panel(
      title = "Raw Data",
      
      layout_sidebar(
        sidebar = sidebar(
          #--- Upload raw data
          tags$h4("Upload Raw Data"),
          tags$h6("When sorted by injection order, the data must begin and end with a QC samples."),
          fileInput(
            inputId = "file1",
            label = "Choose CSV File",
            accept = ".csv",
            buttonLabel = "Browse...",
            placeholder = "No file selected"
          ),
          tags$hr(),
          
          #--- select non-metabolite columns
          tags$h4("Select non-metabolite columns"),
          tags$h6("Please select columns for sample, batch, class, and order."),
          uiOutput("column_selectors"),
          uiOutput("column_warning"),
          #tags$h6("Coming soon: option to select other non-metabolite columns."),
          tags$hr(),
          
          #--- Filter metabolites
          tags$h4("Filter Raw Data"),
          sliderInput(
            inputId = "Frule",
            label = "Acceptable % of missing values per metabolite",
            min = 0,
            max = 100,
            value = 20
          ),
          tags$hr(),
          actionButton(inputId = "next_correction", label = "Choose Correction Settings"),
          
          width = 400,
        ),
        card(
          card_title("Raw Data (first 10 rows)"),
          tableOutput("contents")
        ),
        card(
          card_title("Basic Information"),
          uiOutput("basic_info")
        ),
        card(
          card_title("Filtering Information"),
          uiOutput("filter_info")
        )
      )
    ),
    
    #--- Step 2: Correction settings
    accordion_panel(
      title = "Correction Settings",
      
      layout_sidebar(
        sidebar = sidebar(
          #--- Impute missing values
          tags$h4("Impute Missing Values"),
          radioButtons(
            inputId = "imputeM",
            label = "Imputation method",
            choices = list("metabolite median" = "median", 
                           "metabolite mean" = "mean", 
                           "class-metabolite median" = "class_median", 
                           "class-metabolite mean" = "class_mean",
                           "KNN" = "KNN", 
                           "minimum value" = "min", 
                           "half minimum value" = "minHalf",
                           "zero" = "zero"),
            selected = "median"
          ),
          tags$hr(),
          
          #--- Choose Correction method
          tags$h4("Correction Method"),
          radioButtons(
            inputId = "corMethod",
            label = "Method",
            choices = list("Random Forest Signal Correction" = "QCRFSC", 
                           "Local Polynomial Fit" = "QCRLSC" 
            ),
            selected = "QCRFSC"
          ),
          tags$hr(),
    
          # After correction filtering
          tags$h4("Post-Correction Filtering"),
          sliderInput(
            inputId = "rsd_filter",
            label = "Metabolite RSD% threshold for QC samples",
            min = 0,
            max = 100,
            value = 50
          ),
          tags$hr(),
          
          # After correction scaling / normalization
          tags$h4("Post-Correction Transformation or Normalization"),
          radioButtons(
            inputId = "transform",
            label = "Method",
            choices = list("Log 2 transformation" = "log2", 
                           "Total Ratio Normalization" = "TRN", 
                           "None" = "none"),
            selected = "none"
          ),
          tags$hr(),
          
          actionButton(inputId = "correct", label = "Correct Data"),
          width = 400,
        ),
        card(
          card_title("Correction Information"),
          uiOutput("correction_info") %>% withSpinner(color = "#404040"),
        ),
        card(
          card_title("Corrected Data"),
          tableOutput("cor_data"),
          uiOutput("download_corr_btn", container = div, 
                   style = "position: absolute; bottom: 15px; right: 15px;")
        )
      )
    ),
    
    #--- Step 4: plots
    accordion_panel(
      title = "Evaluation Metrics and Visualization",
      
      layout_sidebar(
        sidebar = sidebar(
          #--- plot metabolite
          tags$h4("Metabolite Scatter plot"),
          uiOutput("met_plot_selectors")
        ),
        width = 400
      ),
      card(
        card_title("Metabolite Scatter Plot"),
        plotOutput("metab_scatter", height = "600px")
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  data_raw <- reactive({
    req(input$file1)
    df <- read.csv(input$file1$datapath, header=TRUE)
  })
  
  #–– 2) preview
  output$contents <- renderTable(head(data_raw(),10))
  
  #–– 3) column selection & warning
  output$column_selectors <- renderUI({
    req(data_raw())
    cols <- names(data_raw())
    
    dropdown_choices <- c("Select a column..." = "", cols)
    
    tagList(
      selectInput("sample_col", "sample column", choices = dropdown_choices, selected = ""),
      selectInput("batch_col", "batch column", choices = dropdown_choices, selected = ""),
      selectInput("class_col", "class column", choices = dropdown_choices, selected = ""),
      selectInput("order_col", "order column", choices = dropdown_choices, selected = "")
    )
  })
  
  output$column_warning <- renderUI({
    req(data_raw())
    selected <- c(
      input$sample_col,
      input$batch_col,
      input$class_col,
      input$order_col
    )
    columnWarningUI(data_raw(), selected)
  })
  
  #–– 4) cleaned data + track replacements
  cleaned <- reactive({
    sel <- c(
      input$sample_col,
      input$batch_col,
      input$class_col,
      input$order_col
    )
    req(!any(sel==""), length(unique(sel))==4)
    cleanData(data_raw(), sel[1], sel[2], sel[3], sel[4])
  })
  
  #–– 5) basic info
  output$basic_info <- renderUI({
    cleaned_data <- cleaned()
    req(cleaned_data)
    basicInfoUI(cleaned_data$df, cleaned_data$replacement_counts)
  })
  
  #–– 6) filter step
  filtered <- reactive({
    cleaned_data <- cleaned()
    req(cleaned_data)
    filter_data(cleaned_data$df, setdiff(names(cleaned_data$df),
                                      c("sample","batch","class","order")), input$Frule)
  })
  output$filter_info <- renderUI({
    filtered_data <- filtered()
    req(filtered_data)
    filterInfoUI(filtered_data$removed_cols)
  })
  
  observeEvent(input$next_correction, {
    accordion_panel_close(id = "main_steps", value = "Raw Data" , session = session)
    accordion_panel_open(id = "main_steps", value = "Correction Settings", session = session)
  })
  
  imputed <- eventReactive(input$correct, {
    filtered_result <- filtered()
    req(filtered_result)
    
    impute_missing(filtered_result$df_filtered,
                   setdiff(names(filtered_result$df_filtered), c("sample", "batch", "class", "order")), 
                   input$imputeM)
  })
  corrected <- eventReactive(input$correct, {
    imputed_result <- imputed()
    req(imputed_result)

    correct_data(imputed_result$df_imputed, 
                 setdiff(names(imputed_result$df_imputed), c("sample", "batch", "class", "order")),
                 input$corMethod)
  })
  
  filtered_corrected <- eventReactive(input$correct, {
    corrected_result <- corrected()
    req(corrected_result)
    
    rsd_filter(corrected_result$df_corrected, input$rsd_filter, c("sample", "batch", "class", "order"))
  })
  
  
  output$correction_info <- renderUI({
    imputed_result <- imputed()
    corrected_result <- corrected()
    filtered_corrected_result <- filtered_corrected()
    req(imputed_result, corrected_result, filtered_corrected_result)
    correctionInfoUI(imputed_result, corrected_result, filtered_corrected_result, input$imputeM, input$corMethod)
    
  })
  
  output$cor_data <- renderTable({
    fil_cor_result <- filtered_corrected()
    req(fil_cor_result)
    head(fil_cor_result$filtered_df, n = 10)
  })
  
  output$download_corr_btn <- renderUI({
    # only run once correction has happened
    req(filtered_corrected())
    
    # show the button
    downloadButton(
      outputId = "download_corr_data",
      label    = "Download CSV",
      class    = "btn btn-primary"
    )
  })
  
  output$download_corr_data <- downloadHandler(
    filename = function() {
      paste0("corrected_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_corrected()$filtered_df, file, row.names = FALSE)
    }
  )
  
  output$met_plot_selectors <- renderUI({
    req(filtered(), filtered_corrected())
    raw_cols <- setdiff(names(filtered()$df_filtered), c("sample", "batch", "class", "order"))
    cor_cols <- setdiff(names(filtered_corrected()$filtered_df), c("sample", "batch", "class", "order"))
    cols <- intersect(raw_cols, cor_cols)
    
    dropdown_choices <- c("Select a column..." = "", cols)
    
    tagList(
      selectInput("met_col", "Metabolite column", choices = dropdown_choices, selected = "")
    )
  })
  output$metab_scatter <- renderPlot({
    req(input$met_col, filtered(), filtered_corrected())
    met_scatter_rf(
      data_raw = filtered()$df_filtered,
      data_cor = filtered_corrected()$filtered_df,
      i = input$metab
    )
  })
  
}


shinyApp(ui = ui, server = server)
# Load packages

library(shiny)
library(bslib)
library(dplyr)
source("helpers.R")


# Define UI

ui <- page_sidebar(
  theme = bs_theme(preset = "cosmo"),
  sidebar = sidebar(
    
    #--- Upload Raw data
    tags$h4("Upload Raw Data"),
    tags$h6("When sorted in order, the data must begin and end with a QC sample."),
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
    tags$h6("Must select unique columns for sample, batch, class, and order."),
    uiOutput("column_selectors"),
    uiOutput("column_warning"),
    tags$hr(),
    
    #--- Filter metabolites
    tags$h4("Filter Raw Data"),
    sliderInput(
      inputId = "missingRule",
      label = "Acceptable % of missing values per metabolite",
      min = 0,
      max = 100,
      value = 20
    ),
    tags$hr(),
    
    #--- impute missing values
    tags$h4("Impute Missing Values"),
    radioButtons(
      inputId = "imputeM",
      label = "Imputation method",
      choices = list("metabolite median" = "median", "metabolite mean" = "met_mean", 
                     "class-metabolite median" = "class_med", "class-metabolite mean" = "class_mean",
                     "KNN" = "KNN", "minimum value" = "min", "half minimum value" = "minHalf",
                     "zero" = "zero"),
      selected = "median"
    ),
    tags$hr(),
    
    #--- Choose Correction method
    tags$h4("Correction Method"),
    radioButtons(
      inputId = "MLmethod",
      label = "Method",
      choices = list("Random Forest Signal Correction" = "QCRFSC", "Local Polynomial Fit" = "QCRLSC" 
                     ),
      selected = "QCRFSC"
    ),
    tags$hr(),
    
    # After correction filtering
    tags$h4("Post-Correction Filtering"),
    sliderInput(
      inputId = "coCV",
      label = "Acceptable metabolite %RSD for QC samples",
      min = 0,
      max = 100,
      value = 50
    ),
    tags$hr(),
    
    width = 400,
  ),
  
  card(
    # Display raw data table
    card_title("Raw Data (first 10 rows)"),
    tableOutput("contents")
  ),
  card(
    card_title("Basic Information"),
    uiOutput("basic_info")
  ),
  card(
    card_title("Post-Corrected Data")
  )
)

# Define server
server <- function(input, output, session) {
  # upload data
  uploaded_data <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = TRUE)
  })
  # display raw data table 
  output$contents <- renderTable({
    head(uploaded_data(), n = 10)
  })
  
  # Get column names and select non-metabolite columns
  output$column_selectors <- renderUI({
    req(uploaded_data())
    cols <- names(uploaded_data())
    
    dropdown_choices <- c("Select a column..." = "", cols)
    
    tagList(
      selectInput("sample_col", "sample column", choices = dropdown_choices, selected = ""),
      selectInput("batch_col", "batch column", choices = dropdown_choices, selected = ""),
      selectInput("class_col", "class column", choices = dropdown_choices, selected = ""),
      selectInput("order_col", "order column", choices = dropdown_choices, selected = "")
    )
  })
  
  # Make sure all are selected and are unique
  output$column_warning <- renderUI({
    selected <- c(
      input$sample_col,
      input$batch_col,
      input$class_col,
      input$order_col
    )
    
    if (any(selected == "")) {
      tags$span(style = "color: darkorange;", "Please select all four columns.")
    } else if (length(unique(selected)) < 4 & length(selected) > 0) {
      tags$span(style = "color: red;", "Each selected column must be unique.")
    } else {
      NULL
    }
  })
  
  replacement_summary <- reactiveVal(NULL)
  # Rename columns, QCs and replace any non-numeric values in metabolite columns with NA
  cleaned_data <- reactive({
    req(uploaded_data())
    req(input$sample_col, input$batch_col, input$class_col, input$order_col)
    
    selected_cols <- c(input$sample_col, input$batch_col, input$class_col, input$order_col)
    if (any(selected_cols == "") || length(unique(selected_cols)) < 4) {
      return(NULL)
    }
    
    df <- uploaded_data()
    colnames(df)[colnames(df) == input$sample_col] <- "sample"
    colnames(df)[colnames(df) == input$batch_col] <- "batch"
    colnames(df)[colnames(df) == input$class_col] <- "class"
    colnames(df)[colnames(df) == input$order_col] <- "order"
    
    # Identify metabolite columns
    metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
    
    # Track replacements
    replacement_counts <- data.frame(
      metabolite = metab_cols,
      replaced_nan_string = integer(length(metab_cols)),
      replaced_zero = integer(length(metab_cols)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(metab_cols)) {
      col <- metab_cols[i]
      original <- df[[col]]
      numeric_version <- suppressWarnings(as.numeric(as.character(original)))
      
      non_numeric_count <- sum(is.na(numeric_version) & !is.na(original))
      zero_count <- sum(numeric_version == 0, na.rm = TRUE)
      
      # Replace values
      numeric_version[numeric_version == 0] <- NA
      df[[col]] <- numeric_version
      
      replacement_counts$non_numeric_replaced[i] <- non_numeric_count
      replacement_counts$zero_replaced[i] <- zero_count
    }
    
    replacement_summary(replacement_counts)
    
    # Normalize QC class labels
    df$class[df$class %in% c("QC", "qc", "Qc")] <- NA
    
    return(df)
  })
  
  output$basic_info <- renderUI({
    req(cleaned_data())
    
    # Get info from data
    df <- cleaned_data()
    metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
    n_metab <- length(metab_cols)
    n_missv <- sum(is.na(df[, metab_cols]))
    n_qcs <- sum(is.na(df$class))
    n_samp <- sum(!is.na(df$class))
    n_bat <- length(unique(df$batch))
    n_class <- length(unique(df$class[!is.na(df$class)]))
    class_list <- sort(unique(df$class[!is.na(df$class)]))
    
    qc_per_batch <- df %>%
      group_by(batch) %>%
      summarise(qc_in_class = sum(is.na(class)), .groups = "drop")
    
    replace_df <- replacement_summary()
    total_replaced <- sum(replace_df$replaced_nan_string + replace_df$replaced_zero)
    
    # Helper for pretty metric cards
    metric_card <- function(label, value, color = "#f8f9fa") {
      tags$div(
        style = paste(
          "background-color:", color,
          "; padding: 10px 15px; border-radius: 8px; margin: 5px;",
          "min-width: 200px; max-width: 220px; flex: 1;"
        ),
        tags$h5(style = "margin: 0;", label),
        tags$p(style = "font-size: 1.4em; font-weight: bold; margin: 0;", value)
      )
    }
    
    # Table of QC samples per batch
    na_table_html <- tags$div(
      style = "margin-top: 15px;",
      tags$h5("Number of QC Samples per Batch"),
      tags$table(
        class = "table table-bordered table-sm",
        tags$thead(
          tags$tr(
            tags$th("Batch"),
            tags$th("QCs in Batch")
          )
        ),
        tags$tbody(
          lapply(1:nrow(qc_per_batch), function(i) {
            tags$tr(
              tags$td(as.character(qc_per_batch$batch[i])),
              tags$td(qc_per_batch$qc_in_class[i])
            )
          })
        )
      )
    )
    
    # Class badges
    class_badges <- tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;",
      lapply(class_list, function(cls) {
        tags$span(
          style = "background-color: #e9ecef; padding: 5px 10px; border-radius: 12px;",
          as.character(cls)
        )
      })
    )
    
    # Main layout
    tagList(
      if (total_replaced > 0) {
        tags$span(style = "color: darkred; font-weight: bold;", 
                  paste(total_replaced, "non-numeric or zero metabolite values were removed."))
      },
      tags$div(
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px;",
        metric_card("Metabolite Columns", n_metab),
        metric_card("Missing Values (Metabolites)", n_missv),
        metric_card("QC Samples", n_qcs),
        metric_card("Samples", n_samp),
        metric_card("Batches", n_bat),
        metric_card("Classes", n_class),
      ),
      tags$div(
        style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 20px;",
        
        # Class badge column
        tags$div(
          style = "flex: 1; min-width: 250px;",
          tags$h5("Unique Classes"),
          class_badges
        ),
        
        # QC per batch table column
        tags$div(
          style = "flex: 1; min-width: 250px;",
          tags$h5("Number of QC Samples per Batch"),
          tags$table(
            class = "table table-bordered table-sm",
            tags$thead(
              tags$tr(
                tags$th("Batch"),
                tags$th("QCs in Batch")
              )
            ),
            tags$tbody(
              lapply(1:nrow(qc_per_batch), function(i) {
                tags$tr(
                  tags$td(as.character(qc_per_batch$batch[i])),
                  tags$td(qc_per_batch$qc_in_class[i])
                )
              })
            )
          )
        )
      )
    )
  })
  
}


shinyApp(ui = ui, server = server)
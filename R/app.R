# Load packages

library(shiny)
library(bslib)
library(dplyr)
source("helpers.R")


# Define UI

ui <- page_fillable(
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
          tags$h6("Please select unique columns for sample, batch, class, and order."),
          uiOutput("column_selectors"),
          uiOutput("column_warning"),
          tags$h6("Coming soon: option to select other non-metabolite columns."),
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
    
    #--- Step 2: Corerction settings
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
            inputId = "coCV",
            label = "Acceptable metabolite %RSD for QC samples",
            min = 0,
            max = 100,
            value = 50
          ),
          tags$hr(),
          actionButton(inputId = "correct", label = "Correct Data"),
          width = 400,
        ),
        card(
          card_title("Correction Results"),
          uiOutput("imputed_info"),
          uiOutput("correction_info"),
          tableOutput("cor_data")
        )
      )
    )
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
        tags$p(style = "font-size: 1.4em; font-weight: bold; margin: 0;", value),
        tags$h5(style = "margin: 0;", label),
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
        style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 10px;",
        # left side
        tags$div(
          style = "display: grid; grid-template-columns: repeat(1, 1fr); gap: 20px;",
          
          # top left
          tags$div(
            style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-top 15px;",
            metric_card("Metabolite Columns", n_metab),
            metric_card("Missing Values", n_missv),
            metric_card("QC Samples", n_qcs),
            metric_card("Samples", n_samp),
            metric_card("Batches", n_bat),
            metric_card("Classes", n_class),
          ),
          
          #bottom left
          tags$div(
            style = "flex: 1; min-width: 250px;",
            tags$h5("Unique Classes"),
            class_badges
          )
          
        ),
        
        # right half
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
  
  filtered_result <- reactive({
    df <- cleaned_data()
    req(df)
    
    metabolite_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
    filter_data(df, metabolite_cols, input$Frule)
  })
  
  output$filter_info <- renderUI({
    result <- filtered_result()
    req(result)
    
    n_removed <- length(result$removed_cols)
    removed_names <- result$removed_cols
    
    if (n_removed == 0) {
      return(tags$span(style = "color: darkgreen; font-weight: bold;",
                       "No metabolite columns were removed based on missing value threshold."))
    }
    
    tagList(
      tags$span(style = "color: darkorange; font-weight: bold;",
                paste(n_removed, "metabolite columns were removed based on the missing value threshold.")),
      tags$br(),
      tags$h5("Removed Columns:"),
      tags$ul(
        lapply(removed_names, function(name) {
          tags$li(name)
        })
      )
    )
  })
  
  observeEvent(input$next_correction, {
    accordion_panel_close(id = "main_steps", value = "Raw Data" , session = session)
    accordion_panel_open(id = "main_steps", value = "Correction Settings", session = session)
  })
  
  imputed_result <- eventReactive(input$correct, {
    result <- filtered_result()
    req(result)
    
    filtered_df <- result$df_filtered
    metab_cols <- setdiff(names(filtered_df), c("sample", "batch", "class", "order"))
    impute_missing(filtered_df, metab_cols, input$imputeM, input$class_col)
  })
  
  output$imputed_info <- renderUI({
    result <- imputed_result()
    req(result)
    
    tags$span(paste(result$n_missv, "missing values imputed with", result$impute_str))
  })
  
  corrected_result <- eventReactive(input$correct, {
    result <- imputed_result()
    req(result)
    
    imputed_df <- result$df_imputed
    metab_cols <- setdiff(names(imputed_df), c("sample", "batch", "class", "order"))
    correct_data(imputed_df, metab_cols, input$corMethod)
  })
  output$correction_info <- renderUI({
    result <- corrected_result()
    req(result)
    tags$span(paste("Data corrected with", result$cor_str))
  })
  output$cor_data <- renderTable({
    result <- corrected_result()
    req(result)
    head(result$df_corrected, n = 10)
  })
}


shinyApp(ui = ui, server = server)
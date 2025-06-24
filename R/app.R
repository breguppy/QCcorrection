# Load packages

library(shiny)
library(bslib)
library(dplyr)


# Define UI

ui <- page_sidebar(
  theme = bs_theme(preset = "cosmo"),
  sidebar = sidebar(
    
    # Upload Raw data
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
    
    # select non-metabolite columns
    tags$h4("Select non-metabolite columns"),
    tags$h6("Must select unique columns for sample, batch, class, and order."),
    uiOutput("column_selectors"),
    uiOutput("column_warning"),
    tags$hr(),
    
    # Filter metabolites
    tags$h4("Filter Raw Data"),
    sliderInput(
      inputId = "missingRule",
      label = "Acceptable % of missing values per metabolite",
      min = 0,
      max = 100,
      value = 20
    ),
    tags$hr(),
    
    # impute missing values
    tags$h4("Impute Missing Values"),
    radioButtons(
      inputId = "imputeM",
      label = "Imputation method",
      choices = list("metabolite median" = 1, "metabolite mean" = 2, 
                     "class-metabolite median" = 3, "class-metabolite mean" = 4,
                     "KNN" = 5, "minimum value" = 6, "half minimum value" = 7,
                     "zero" = 8),
      selected = 1
    ),
    tags$hr(),
    
    # impute missing values
    tags$h4("Correction Method"),
    radioButtons(
      inputId = "MLmethod",
      label = "Method",
      choices = list("Random Forest Signal Correction" = 1, "Local Polynomial Fit" = 2 
                     ),
      selected = 1
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
    card_title("Raw Data (first 10 rows)"),
    tableOutput("contents")
  ),
  card(
    card_title("Basic Information"),
    uiOutput("basic_info")
  )
)

# Define server
server <- function(input, output, session) {
  # Display raw data table
  uploaded_data <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, header = TRUE)
  })
  
  output$contents <- renderTable({
    head(uploaded_data(), n = 10)
  })
  
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
      x <- df[[col]]
      
      # Count "NaN" strings
      replaced_nan <- sum(x == "NaN", na.rm = TRUE)
      
      # Convert to numeric (handles characters and factors)
      x <- as.numeric(as.character(x))
      
      # Count exact 0s (after conversion)
      replaced_zero <- sum(x == 0, na.rm = TRUE)
      
      # Replace
      x[x == 0] <- NA
      df[[col]] <- x
      df[[col]][df[[col]] == "NaN"] <- NA  # just in case some slipped through
      
      # Save replacement info
      replacement_counts$replaced_nan_string[i] <- replaced_nan
      replacement_counts$replaced_zero[i] <- replaced_zero
    }
    
    # Save summary table to global reactive
    replacement_summary(replacement_counts)
    
    
    df$class[df$class == "QC"] <- NA
    df$class[df$class == "qc"] <- NA
    df$class[df$class == "Qc"] <- NA    
    return(df)
  })
  
  output$basic_info <- renderUI({
    req(cleaned_data())
    
    df <- cleaned_data()
    metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
    n_metab <- length(metab_cols)
    n_missv <- sum(is.na(df[, metab_cols]))
    n_qcs <- sum(is.na(df$class))
    n_samp <- sum(!is.na(df$class))
    n_bat <- length(unique(df$batch))
    n_class <- length(unique(df$class[!is.na(df$class)]))
    class_list <- sort(unique(df$class[!is.na(df$class)]))
    
    # Compute NA count in 'class' per batch
    qc_per_batch <- df %>%
      group_by(batch) %>%
      summarise(qc_in_class = sum(is.na(class)), .groups = "drop")
    
    # Convert the table to HTML
    na_table_html <- tagList(
      tags$strong("Number of QCs per batch:"),
      tags$table(
        class = "table table-sm",
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
    
    replace_df <- replacement_summary()
    total_replaced <- sum(replace_df$replaced_nan_string + replace_df$replaced_zero)
    
    tagList(
      tags$p(paste("Number of metabolite columns:", n_metab)),
      tags$p(paste("Number of QC samples:", n_qcs)),
      tags$p(paste("Number of samples:", n_samp)),
      tags$p(paste("Number of batches:", n_bat)),
      tags$p(paste("Number of classes:", n_class)),
      tags$p("Classes:"),
      tags$div(
        style = "display: flex; flex-wrap: wrap; gap: 12px;",
        lapply(class_list, function(cls) {
          tags$div(style = "width: 25%;", as.character(cls))
        })
      ),
      tags$p(paste("Number of missing values:", n_missv)),
      tags$p(paste("Total values changed in metabolite columns:", total_replaced)),
      na_table_html
    )
  })
  
}


shinyApp(ui = ui, server = server)
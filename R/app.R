# Load packages

library(shiny)
library(bslib)


# Define UI

ui <- page_sidebar(
  theme = bs_theme(preset = "cosmo"),
  sidebar = sidebar(
    
    # Upload Raw data
    tags$h3("Upload Raw Data"),
    tags$h6("When sorted in order, the data must begin and end with a QC sample."),
    fileInput(
      inputId = "file1",
      label = "Choose CSV File",
      accept = ".csv",
      buttonLabel = "Browse...",
      placeholder = "No file selected",
    ),
    tags$hr(),
    
    # select non-metabolite columns
    tags$h3("Select non-metabolite columns"),
    tags$h6("Must select unique columns for sample, batch, class, and order."),
    uiOutput("column_selectors"),
    uiOutput("column_warning"),
    tags$hr(),
    
    # Filter metabolites
    tags$h3("Filter Data"),
    sliderInput(
      inputId = "missingRule",
      label = "Acceptable % of missing values per metabolite",
      min = 0,
      max = 100,
      value = 20
    ),
    tags$hr(),
    
    # After correction filtering
    tags$h3("Post-Correction Filtering"),
    sliderInput(
      inputId = "coCV",
      label = "Acceptable metabolite %RSD for QC samples",
      min = 0,
      max = 100,
      value = 50
    ),
    tags$hr(),
    
    width = 350,
  ),
  
  card(
    tableOutput("contents")
  )
)
# Define server

server <- function(input, output, session) {
  
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
  #input$sample_col  # selected Sample column
  #input$batch_col   # selected Batch column
  #input$class_col   # selected Class column
  #input$order_col 
}
# Create a Shiny app object

shinyApp(ui = ui, server = server)
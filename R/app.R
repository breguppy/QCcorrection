# Load packages

library(shiny)
library(bslib)


# Define UI

ui <- page_sidebar(
  theme = bs_theme(preset = "minty"),
  sidebar = sidebar(
      fileInput(
        inputId = "file1",
        label = "Choose CSV File",
        accept = ".csv",
        width = NULL,
        buttonLabel = "Browse...",
        placeholder = "No file selected",
        capture = NULL
      )
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
}
# Create a Shiny app object

shinyApp(ui = ui, server = server)
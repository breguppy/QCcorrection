# mod_export.R

mod_export_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "4. Export Corrected Data, Plots, and Report", 
    value = "tab_export",
    card(
    card_title("Download Data, Plots, and Report"),
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
             )),
      column(4, tags$span(icon("file-pdf"), " correction_report.pdf"),
             tags$ul(
               tags$li("Report describing the correction steps and figures for evaluating the correction process.")
             ))
    ),
    uiOutput(ns("download_all_ui"))
  ))
}

mod_export_server <- function(id, filtered, filtered_corrected, transformed, params) {
  moduleServer(id, function(input, output, session) {
  })
}
#' Export module
#' @keywords internal
#' @noRd

mod_export_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "4. Export Corrected Data, Plots, and Report", 
    value = "tab_export",
    card(
    card_title("Download Data, Plots, and Report"),
    tags$span("Download all to get a ", icon("folder"), " zipped folder containing:"),
    fluidRow(
      column(3, tags$span(icon("file-excel"), "corrected_data_*today's_date*.xlsx"),
             tags$ul(
               tags$li("0. Raw Data"),
               tags$li("1. Correction Settings"),
               tags$li("2. Drift Normalized"),
               tags$li("3. Scaled or Normalized"),
               tags$li("4. Grouped Data Organized"),
               tags$li("5. Grouped Data Fold Change (only when provided a control class)"),
               tags$li("Appendix1. Metaboanalyst Ready")
             )),
      column(2, tags$span(icon("file-excel"), "rsd_stats_*today's_date*.xlsx"),
             tags$ul(
               tags$li("Raw RSD"),
               tags$li("Corrected RSD or Transformed Corrected RSD"),
               tags$li("RSD Comparison")
               )),
      column(2, tags$span(icon("file-excel"), "candidate_outliers_*today's_date*.xlsx"),
             tags$ul(
               tags$li("QC RSD"),
               tags$li("Sample MD"),
               tags$li("Candidates"),
               tags$li("Confirmations")
               )),
      column(2, tags$span(icon("folder"), " figures"),
             tags$ul(
               tags$li(icon("folder"), " metabolite figures"),
               tags$li(icon("folder"), " RSD figures"),
               tags$li(icon("folder"), "PCA plots")
             )),
      column(3, tags$span(icon("file-pdf"), " correction_report.pdf"),
             tags$ul(
               tags$li("Report describing the correction steps and figures for evaluating the correction process.")
             ))
    ),
    uiOutput(ns("download_all_ui"))
  ))
}

mod_export_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    d <- reactive(data())
    p <- reactive(params())
    
    output$download_all_ui <- renderUI({
      req(d()$transformed)
      downloadButton(
        outputId = ns("download_all_zip"),
        label = "Download All",
        class = "btn-primary btn-lg"
      )
    })
    
    #-- Allow user to download corrected data, figures, and correction report.
    output$download_all_zip <- downloadHandler(
      filename = function() {
        sprintf("corrected_data_plots_report_%s.zip", Sys.Date())
      },
      content = function(file) {
        .require_pkg("zip", "create a zip archive")
        base_dir <- tempfile("bundle_")
        dir.create(base_dir)
        on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)
        
        # create and save corrected data file
        xlsx_path <- file.path(base_dir, sprintf("corrected_data_%s.xlsx", Sys.Date()))
        wb <- export_xlsx(p(), d())
        saveWorkbook(wb, xlsx_path, overwrite = TRUE)
        
        # Create and save rsd stats data file
        stats_xlsx_path <- file.path(base_dir, sprintf("rsd_stats_%s.xlsx", Sys.Date()))
        stats_wb <- export_stats_xlsx(p(), d())
        saveWorkbook(stats_wb, stats_xlsx_path, overwrite = TRUE)
        
        # Create and save outlier data file
        outlier_xlsx_path <- file.path(base_dir, sprintf("candidate_outliers_%s.xlsx", Sys.Date()))
        outlier_wb <- export_outliers_xlsx(p(), d())
        saveWorkbook(outlier_wb, outlier_xlsx_path, overwrite = TRUE)
        
        # create and save figure folder
        figs <- export_figures(p(), d(), out_dir = base_dir)
        
        # make pdf report
        render_report(p(), d(), out_dir = base_dir)
        
        # make zip file
        rel <- c(
          "figures",                                  
          basename(xlsx_path),
          basename(stats_xlsx_path),
          basename(outlier_xlsx_path),
          "correction_report.html", "correction_report.pdf" 
        )
        rel <- rel[file.exists(file.path(base_dir, rel))]
        
        tmpzip <- tempfile(fileext = ".zip")
        on.exit(unlink(tmpzip, force = TRUE), add = TRUE)
        zip::zipr(zipfile = tmpzip, files = rel, root = base_dir)
        
        file.copy(tmpzip, file, overwrite = TRUE)
      }
    )
    
  })
}
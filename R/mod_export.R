#' @keywords internal

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
        paste0("corrected_data_plots_report_", Sys.Date(), ".zip")
      },
      content = function(file) {
        base_dir <- tempfile("bundle_")
        dir.create(base_dir)
        
        # create and save corrected data file
        cor_data_filename <- paste0("corrected_data_", Sys.Date(), ".xlsx")
        cor_data_path <- file.path(base_dir, cor_data_filename)
        wb <- corrected_file_download(p(), d())
        saveWorkbook(wb, cor_data_path, overwrite = TRUE)
        
        # 2. create and save figure folder
        fig_info <- figure_folder_download(p(), d())
        figs_src <- fig_info$fig_dir
        figs_dir <- file.path(base_dir, "figures"); dir.create(figs_dir)
        
        files <- list.files(figs_src, recursive = TRUE, full.names = TRUE)
        rel   <- list.files(figs_src, recursive = TRUE, full.names = FALSE)
        
        # build target paths
        targets <- file.path(figs_dir, rel)
        
        # create subdirectories (skip ".")
        dirs <- unique(dirname(targets))
        dirs <- dirs[dirs != "."]   
        for (dir in dirs) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        
        # copy files into same structure
        file.copy(from = files, to = targets, overwrite = TRUE)
        
        # 3. make pdf report
        generate_cor_report(p(), d(), base_dir)
        
        # make zip file
        zipfile <- tempfile(fileext = ".zip")
        on.exit(unlink(c(base_dir, zipfile), recursive = TRUE), add = TRUE)
        zip::zipr(zipfile = zipfile,
                  files = list.files(base_dir, full.names = TRUE), 
                  include_directories = TRUE)
        file.copy(zipfile, file, overwrite = TRUE)
      }
    )
    
  })
}
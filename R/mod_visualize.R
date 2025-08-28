# mod_visualize.R

mod_visualize_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "3. Evaluation Metrics and Visualization",
    value = "tab_visualize",
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("3.1 Visualize Correction with Metabolite Scatter Plots"),
        uiOutput(ns("met_plot_selectors")),
        width = 400,
      ),
      plotOutput(ns("metab_scatter"), height = "600px", width = "600px"),
    )),
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("3.2 RSD Evaluation"),
        tags$h6("Evaluate correction method by the change in relative standard deviation (RSD)."),
        radioButtons(ns("rsd_compare"), "Compare raw data to", list("Corrected data" = "filtered_cor_data", "Transformed and corrected data" = "transformed_cor_data"), "filtered_cor_data"),
        radioButtons(ns("rsd_cal"), "Calculate RSD by", list("Metabolite" = "met", "Class and Metabolite" = "class_met"), "met"),
        width = 400,
      ),
      plotOutput(ns("rsd_comparison_plot"), height = "540px", width = "900px")
    )),
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("3.3 PCA Evaluation"),
        tags$h6("Evaluate correction using principal component analysis (PCA)."),
        radioButtons(ns("pca_compare"), "Compare raw data to", list("Corrected data" = "filtered_cor_data", "Transformed and corrected data" = "transformed_cor_data"), "filtered_cor_data"),
        radioButtons(ns("color_col"), "Color PCA by", list("batch" = "batch", "class" = "class"), "batch"),
        width = 400,
      ),
      plotOutput(ns("pca_plot"), height = "530px", width = "1000px")
    )),
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("3.4 Download Figures Only"),
        tooltip(
          radioButtons(ns("fig_format"), "Select figure format:", c("PDF" = "pdf", "PNG" = "png"), "pdf"),
          "All figures will be saved in this format after clicking download button here or on tab 4. Export Corrected Data, Plots, and Report", "right"
        ),
        width = 400,
      ),
      uiOutput(ns("download_fig_zip_btn")),
      uiOutput(ns("progress_ui")),
    )),
    card(
      actionButton(ns("next_export"), "Next: Export Corrected Data and Plots", class = "btn-primary btn-lg"),
    )
  )
}

mod_visualize_server <- function(id, filtered, filtered_corrected, transformed, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$rsd_comparison_plot <- renderPlot(execOnResize = FALSE, res = 120,{
      req(input$rsd_compare, input$rsd_cal)
      
      make_rsd_plot(input, list(filtered = filtered(), filtered_correced = filtered_corrected(), transformed = transformed()))
    })
    
    #-- PCA plot
    output$pca_plot <- renderPlot({
      req(input$pca_compare, input$color_col)
      make_pca_plot(input, rv)
    }, res = 120)
    
    #-- Let user select which metabolite to display in scatter plot
    output$met_plot_selectors <- renderUI({
      req(filtered(), filtered_corrected())
      raw_cols <- setdiff(names(filtered()$df),    c("sample","batch","class","order"))
      cor_cols <- setdiff(names(filtered_corrected()$df), c("sample","batch","class","order"))
      cols <- intersect(raw_cols, cor_cols)
      validate(need(length(cols) >= 1, "No overlapping metabolites between raw and corrected data."))
      selectInput(ns("met_col"), "Metabolite column", choices = cols, selected = cols[1])
    })
    
    output$metab_scatter <- renderPlot({
      req(input$met_col)
      make_met_scatter(rv, input$met_col)
    }, res = 120)
    
    #-- Download all figures as zip folder.
    output$download_fig_zip_btn <- renderUI({
      req(transformed())
      
      div(
        style = "max-width: 300px; display: inline-block;",
        downloadButton(ns("download_fig_zip"), "Download All Figures (Optional)", class = "btn-primary")
      )
    })
    # -- progress bar
    progress_reactive <- reactiveVal(0)
    #-- progress for downloading all images
    output$progress_ui <- renderUI({
      req(progress_reactive() > 0, progress_reactive() <= 1)
      
      div(
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
        figs <- figure_folder_download(input, rv)
        
        zipfile <- tempfile(fileext = ".zip")
        old_wd <- setwd(figs$tmp_dir)
        on.exit({
          unlink(figs$fig_dir, recursive = TRUE)
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
    
    #-- Move to next tab after inspecting the corrected data figures
    observeEvent(input$next_export, {
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_export")
    })
    
    list(progress = progress_reactive)
  })
}
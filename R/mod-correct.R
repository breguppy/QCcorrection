#' Correction module
#'
#' @keywords internal
#' @noRd

mod_correct_ui <- function(id) { 
  ns <- NS(id); 
  nav_panel(
    title = "2. Correction Settings",
    value = "tab_correct",
    card(
      style = "background-color:#eee;",
      tags$h4("2.1 Choose Correction Settings"),
      uiOutput(ns("qc_missing_value_warning")),
      fluidRow(
        column(3, tags$h5("Impute Missing QC Values"), uiOutput(ns("qcImpute"))),
        column(3, tags$h5("Impute Missing Sample Values"), uiOutput(ns("sampleImpute"))),
        column(3, tags$h5("Choose Correction Method"), uiOutput(ns("correctionMethod"))),
        column(3, tags$h5("Unavailable Options"), uiOutput(ns("unavailable_options"))),
        actionButton(ns("correct"), "Correct Data with Selected Settings",
                     class="btn-primary btn-lg", width="100%"),
        div(style="margin:12px 0 0 0;", withSpinner(uiOutput(ns("cor_spinner")),
                                                    color="#404040", size=0.6, proxy.height="22px"))
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.2 Post-Correction Filtering",
          ui_post_cor_filter(ns),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Corrected RSD Summary",
            uiOutput(ns("download_cor_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("RSD summary before and after correction for samples and QCs. ",
                     "RSD summary can also be downloaded on tab 4. Export All"),
            width = 400,
            position = "right"
          ),
        uiOutput(ns("post_cor_filter_info")) %>% withSpinner(color = "#404040") 
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.3 Post-Correction Transformation",
          ui_post_cor_transform(ns),
          uiOutput(ns("trn_withhold_ui")),
          uiOutput(ns("trn_withhold_selectors_ui")),
          width = 400
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
             title = "Download Transformed RSD Summary",
          uiOutput(ns("download_tc_rsd_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          help = c("RSD summary before and after correction and transformation for samples and QCs. ",
                   "RSD summary can also be downloaded on tab 4. Export All"),
          width = 400,
          position = "right"
          ),
         ui_table_scroll("cor_data", ns) %>% withSpinner(color = "#404040")
        )
      )
    ),
    card(
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.4 Candidate Extreme Values",
          ui_detect_outliers_options(ns),
          help = c("The samples listed in the table are considered extreme values by robust z\u002Dscore with a cutoff weighted by QC variability confired with a test chosen by group size or by Mahalanobis distance in PCA.")
        ),
        layout_sidebar(
          sidebar = ui_sidebar_block(
            title = "Download Extreme Value Summary",
            uiOutput(ns("download_ev_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
            help = c("Extreme value summary can also be downloaded on tab 4. Export All"),
            width = 400,
            position = "right"
          ),
          uiOutput(ns("outliers_table"))
        )
      )
    ),
    card(
      #style = "background-color: #eeeeee;",
      layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "2.5 Identify Control Group",
          tooltip(
            checkboxInput(ns("no_control"), "No control group", FALSE),
            "Check the box if The data does not have a control group.", 
            placement = "right"
          ),
          conditionalPanel(
          condition = sprintf("!input['%s']", ns("no_control")),
          uiOutput(ns("control_class_selector"))
          ),
          width = 400
        ),
        layout_sidebar(
        sidebar = ui_sidebar_block(
          title = "Download Corrected and Transformed Data",
          tooltip(
            checkboxInput(ns("keep_corrected_qcs"), "Include QCs in corrected data file", FALSE),
            "Check the box if you want corrected QC values in the downloaded corrected data file.", 
            placement = "right"
          ),
          uiOutput(ns("download_corr_btn"), container = div, style = "position: absolute; bottom: 15px; right: 15px;"),
          tags$h6("Corrected data can also be downloaded on tab 4. Export All"),
          width = 400,
          position = "right")
        ))
    ),
    card(actionButton(ns("next_visualization"), "Next: Evaluate and Visualize Correction",
                      class="btn-primary btn-lg"))
  )}

mod_correct_server <- function(id, data, params) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    d <- reactive(data()) 
    
    filtered_r <- reactive({
      req(d()$filtered)
      d()$filtered
    })
    cleaned_r <- reactive({
      d()$cleaned
    })
    
    output$qc_missing_value_warning <- renderUI({
      df <- filtered_r()$df
      ui_qc_missing_warning(df)
    })
    output$qcImpute <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_qc_impute(df, mc, ns = session$ns)
    })
    output$sampleImpute <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_sample_impute(df, mc, ns = session$ns)
    })
    output$correctionMethod <- renderUI({
      ui_correction_method(filtered_r()$df, ns = session$ns)
    })
    output$unavailable_options <- renderUI({
      df <- filtered_r()$df
      mc <- setdiff(names(df), c('sample','batch','class','order'))
      ui_unavailable_options(df, mc)
    })
    
    
    metab_cols_r <- reactive({
      setdiff(names(filtered_r()$df), c("sample","batch","class","order"))
    })
    
    has_qc_na_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class == "QC")[, mc, drop = FALSE]))
    })
    has_sam_na_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()
      any(is.na(dplyr::filter(df, .data$class != "QC")[, mc, drop = FALSE]))
    })
    
    imputed_r <- reactive({
      df <- filtered_r()$df; mc <- metab_cols_r()
      qcImpute  <- if (!has_qc_na_r())  "nothing_to_impute" else input$qcImputeM
      samImpute <- if (!has_sam_na_r()) "nothing_to_impute" else input$samImputeM
      impute_missing(df, mc, qcImpute, samImpute)
    })
    
    corrected_r <- eventReactive(input$correct, {
      imputed <- isolate(imputed_r()); mc <- isolate(metab_cols_r())
      correct_data(imputed$df, mc, isolate(input$corMethod))
    })
    
    observeEvent(input$correct, ignoreInit = TRUE, {
      shinyjs::disable("correct")
      output$cor_spinner <- renderUI({
        on.exit(shinyjs::enable("correct"), add = TRUE)
        corrected_r(); NULL
      })
    })
    output$cor_spinner <- renderUI(NULL)
    
    filtered_corrected_r <- reactive({
      req(filtered_r(), corrected_r())
      df_corrected <- corrected_r()$df
      
      if (isTRUE(input$remove_imputed)) {
        fil_cor_df <- remove_imputed_from_corrected(filtered_r()$df, df_corrected)
      } else  {
        fil_cor_df <- df_corrected
      }
      
      if (isTRUE(input$post_cor_filter)) {
        fil_cor_df <- filter_by_qc_rsd(fil_cor_df, Inf, c("sample","batch","class","order"))
      } else {
        fil_cor_df <- filter_by_qc_rsd(fil_cor_df, input$rsd_filter, c("sample","batch","class","order"))
      }
      
      fil_cor_df 
    })
    
    output$download_cor_rsd_btn <- renderUI({
      req(filtered_corrected_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_cor_rsd_data"),
            label    = "Download Corrected RSD Summary",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    output$download_cor_rsd_data <- downloadHandler(
      filename = function() {
        sprintf("corrected_rsd_stats_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        p <- list(
          rsd_compare = "filtered_cor_data"
        )
        
        d <- list(
          filtered_corrected = filtered_corrected_r(),
          filtered           = filtered_r()
        )
        
        stats_wb <- export_stats_xlsx(p, d)
        openxlsx::saveWorkbook(stats_wb, file, overwrite = TRUE)
      }
    )
    
    
    transformed_r <- reactive({
      req(filtered_corrected_r())
      withheld <- character(0)
      if (isTRUE(input$trn_withhold_checkbox) && !is.null(input$trn_withhold_n)) {
        for (i in seq_len(input$trn_withhold_n)) {
          col <- input[[paste0("trn_withhold_col_", i)]]
          if (!is.null(col) && col %in% names(filtered_corrected_r()$df)) withheld <- c(withheld, col)
        }
      }
      transform_data(filtered_corrected_r()$df, input$transform, withheld, input$ex_ISTD)
    })
    
    observe({
      req(corrected_r(), input$trn_withhold_checkbox)
      
      max_withhold <- max(ncol(corrected_r()$df) - 4, 0)
      
      output$trn_withhold_ui <- renderUI({
        if (input$transform == "TRN") {
          numericInput(
            inputId = ns("trn_withhold_n"),
            label = "Number of columns to withold from TRN",
            value = 1,
            min = 1,
            max = max_withhold
          )
        }
      })
    })
    output$trn_withhold_selectors_ui <- renderUI({
      req(corrected_r(), input$trn_withhold_n, input$ex_ISTD)
      cols <- names(corrected_r()$df)
      cols <- setdiff(cols, c("sample", "batch", "class", "order"))
      if (input$ex_ISTD) {
        istd <- grep("ISTD", cols, value = TRUE)
        itsd <- grep("ITSD", cols, value = TRUE)
        cols <- setdiff(cols, c(istd, itsd))
      }
      dropdown_choices <- c("Select a column..." = "", cols)
      
      # Generate list of columns to withhold
      lapply(seq_len(input$trn_withhold_n), function(i) {
        selectInput(
          inputId = ns(paste0("trn_withhold_col_", i)),
          label = paste("Select column to withhold #", i),
          choices = dropdown_choices,
          selected = ""
        )
      })
    })
    
    output$post_cor_filter_info <- renderUI({
      req(filtered_corrected_r())
      ui_postcor_filter_info(filtered_corrected_r(), input$rsd_filter, input$post_cor_filter)
    })
    output$cor_data <- renderTable({
      req(transformed_r())
      transformed_r()$df
    })
    
    output$download_tc_rsd_btn <- renderUI({
      req(transformed_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_tc_rsd_data"),
            label    = "Download Transformed RSD Summary",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    output$download_tc_rsd_data <- downloadHandler(
      filename = function() {
        sprintf("transformed_rsd_stats_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        p <- list(
          rsd_compare = "transformed_cor_data"
        )
        
        d <- list(
          filtered_corrected = filtered_corrected_r(),
          filtered           = filtered_r(),
          transformed        = transformed_r()
        )
        
        stats_wb <- export_stats_xlsx(p, d)
        openxlsx::saveWorkbook(stats_wb, file, overwrite = TRUE)
      }
    )
    
    output$outliers_table <- renderUI({
      req(filtered_corrected_r(), transformed_r())
      d <- list(filtered_corrected = filtered_corrected_r(), transformed = transformed_r())
      p <- list(out_data = input$out_data, sample_grouping = input$sample_grouping)
      ui_outliers(p, d)
    })
    
    # add download button for extreme values.
    output$download_ev_btn <- renderUI({
      req(transformed_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_ev_data"),
            label    = "Download Extreme Value Summary",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    output$download_ev_data <- downloadHandler(
      filename = function() {
        sprintf("extreme_values_%s.xlsx", Sys.Date())
      },
      content = function(file) {
        d <- list(filtered_corrected = filtered_corrected_r(), transformed = transformed_r())
        p <- list(out_data = input$out_data, sample_grouping = input$sample_grouping)
        
        outlier_wb <- export_outliers_xlsx(p, d)
        openxlsx::saveWorkbook(outlier_wb, file, overwrite = TRUE)
      }
    )
    
    output$control_class_selector <- renderUI({
      req(cleaned_r())
      df <- cleaned_r()$df
      classes <- unique(df$class[df$class != "QC"])
      dropdown_choices <- c("Select a class..." = "", classes)
      tooltip(
        selectInput(
          ns("control_class"),
          "Control Group",
          choices = dropdown_choices,
          selected = ""
        ),
        "Name of control samples in class column. This class's average will be used to compute fold changes in the corrected data file.",
        placement = "right"
      )
    })
    
    
    # button for downloading corrected data.
    output$download_corr_btn <- renderUI({
      req(transformed_r())
      
      div(
        style = "width: 100%; text-align: center;",
        div(
          style = "max-width: 250px; display: inline-block;",
          downloadButton(
            outputId = ns("download_corr_data"),
            label    = "Download Corrected and Transformed Data",
            class    = "btn btn-secondary"
          )
        )
      )
    })
    
    output$download_corr_data <- downloadHandler(
      filename = function() {
        paste0("corrected_data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        fc <- isolate(filtered_corrected_r())
        tr <- isolate(transformed_r())
        cr <- isolate(corrected_r())
        p_in <- params()  
        
        p <- list(
          sample_col        = p_in$sample_col,
          batch_col         = p_in$batch_col,
          class_col         = p_in$class_col,
          order_col         = p_in$order_col,
          Frule             = p_in$Frule,
          remove_imputed    = isTRUE(input$remove_imputed),
          rsd_cutoff        = fc$rsd_cutoff,
          transform         = input$transform,
          ex_ISTD           = isTRUE(input$ex_ISTD),
          keep_corrected_qcs= isTRUE(input$keep_corrected_qcs),
          no_control        = isTRUE(input$no_control),
          control_class     = input$control_class %||% ""
        )
        
        rv <- list(
          cleaned            = cleaned_r(),
          filtered           = filtered_r(),
          imputed            = imputed_r(),
          corrected          = cr,
          filtered_corrected = fc,
          transformed        = tr
        )
        
        wb <- export_xlsx(p, rv)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    observeEvent(input$next_visualization, {
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_visualize")
    })
    
    correct_params <- reactive(list(
      qcImputeM          = input$qcImputeM %||% "median",
      samImputeM         = input$samImputeM %||% "median",
      remove_imputed     = isTRUE(input$remove_imputed),
      post_cor_filter    = input$post_cor_filter,
      rsd_cutoff         = filtered_corrected_r()$rsd_cutoff,
      transform          = input$transform,
      ex_ISTD            = isTRUE(input$ex_ISTD),
      out_data           = input$out_data,
      sample_grouping    = input$sample_grouping,
      keep_corrected_qcs = isTRUE(input$keep_corrected_qcs),
      no_control         = isTRUE(input$no_control),
      control_class      = input$control_class %||% ""
    ))
    
    list(
      imputed            = imputed_r,
      corrected          = corrected_r,
      filtered_corrected = filtered_corrected_r,
      transformed        = transformed_r,
      params             = correct_params
    )
  })
}
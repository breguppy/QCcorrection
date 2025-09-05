#' @keywords internal 

mod_import_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "1. Import Raw Data",
    value = "tab_import",
    card(full_screen = TRUE,
         layout_sidebar(
           sidebar = sidebar(
             tags$h4("1.1 Upload Raw Data"),
             fileInput(ns("file1"), "Choose Raw Data File (.csv, .xls, or .xlsx)",
                       accept = c(".csv",".xls",".xlsx"), 
                       buttonLabel = "Browse...",
                       placeholder = "No file selected"),
             tags$h6("Raw data must be on the first sheet of .xls or .xlsx file."),
             tags$h6("Data must begin and end with QC samples when sorted by injection order."),
             width = 400
           ),
           div(style="overflow:auto; max-height:400px;", tableOutput(ns("contents")))
         )
    ),
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("1.2 Select non-metabolite columns"),
        tags$h6("Please select columns for sample, batch, class, and order."),
        uiOutput(ns("column_selectors")),
        uiOutput(ns("column_warning")),
        tooltip(
          checkboxInput(ns("withhold_cols"), "Withhold additional columns from correction?", FALSE),
          "Select if there are extra non-metabolite or specific metabolite columns to withhold from correction.", "right"
        ),
        uiOutput(ns("n_withhold_ui")),
        uiOutput(ns("withhold_selectors_ui")),
        width = 400
      ),
      uiOutput(ns("basic_info"))
    )),
    card(layout_sidebar(
      sidebar = sidebar(
        tags$h4("1.3 Filter Raw Data"),
        tooltip(
          sliderInput(ns("Frule"), "Acceptable % missing per metabolite", 0, 100, 20),
          "Metabolites above this missing % are removed.", "right"
        ),
        width = 400
      ),
      uiOutput(ns("filter_info"))
    )),
    card(actionButton(ns("next_correction"), "Next: Choose Correction Settings",
                      class = "btn-primary btn-lg"))
  )
}

mod_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    data_raw <- reactive({
      req(input$file1)
      read_raw_data(input$file1$datapath)
    })
    output$contents <- renderTable(data_raw())
    
    selections_r <- reactive({
      list(
        sample = input$sample_col %||% "",
        batch  = input$batch_col  %||% "",
        class  = input$class_col  %||% "",
        order  = input$order_col  %||% ""
      )
    }) %>% debounce(200)
    
    output$column_selectors <- renderUI({
      req(data_raw()); ui_nonmet_cols(names(data_raw()), ns = NS(session$ns(NULL)))
    })
    output$column_warning <- renderUI({
      req(data_raw()); sel <- selections_r()
      ui_column_warning(data_raw(), c(sel$sample, sel$batch, sel$class, sel$order))
    })
    
    withheld_ids_r <- reactive({
      if (!isTRUE(input$withhold_cols)) return(character(0))
      n <- input$n_withhold %||% 0
      if (n <= 0) return(character(0))
      paste0("withhold_col_", seq_len(n))
    })
    withheld_r <- reactive({
      ids <- withheld_ids_r(); if (!length(ids)) return(character(0))
      vals <- vapply(ids, function(id) input[[id]] %||% "", character(1))
      vals <- unique(vals[nzchar(vals)])
      sel <- selections_r()
      setdiff(vals, c(sel$sample, sel$batch, sel$class, sel$order))
    }) %>% debounce(200)
    
    observe({
      req(data_raw())
      max_withhold <- max(ncol(data_raw()) - 4, 0)
      output$n_withhold_ui <- renderUI({
        if (isTRUE(input$withhold_cols))
          numericInput(ns("n_withhold"), "Number of columns to withhold", 1, 1, max_withhold)
      })
    })
    
    output$withhold_selectors_ui <- renderUI({
      req(data_raw(), input$n_withhold)
      sel <- selections_r()
      cols <- setdiff(names(data_raw()), c(sel$sample, sel$batch, sel$class, sel$order))
      ids <- withheld_ids_r(); if (!length(ids)) return(NULL)
      prev_all <- isolate(vapply(ids, function(id) input[[id]] %||% "", character(1)))
      lapply(seq_along(ids), function(i) {
        id <- ids[i]; prev <- prev_all[i]
        other <- setdiff(prev_all, prev)
        choices_i <- c("Select a column..." = "", setdiff(cols, other))
        selectInput(ns(id), paste("Select column to withhold #", i),
                    choices = choices_i, selected = if (nzchar(prev) && prev %in% choices_i) prev else "")
      })
    })
    
    cleaned_r <- reactive({
      df  <- req(data_raw())
      sel <- selections_r()
      withheld <- withheld_r()
      req(all(nzchar(c(sel$sample, sel$batch, sel$class, sel$order))))
      req(length(unique(c(sel$sample, sel$batch, sel$class, sel$order))) == 4)
      clean_data(df, sel$sample, sel$batch, sel$class, sel$order, withheld)
    }) %>% bindCache(reactiveVal(NULL)(), selections_r(), withheld_r())
    
    output$basic_info <- renderUI({
      cd <- cleaned_r(); req(cd)
      ui_basic_info(cd$df, cd$replacement_counts)
    })
    
    filtered_r <- reactive({
      cd <- req(cleaned_r())
      filter_by_missing(cd$df, setdiff(names(cd$df), c("sample","batch","class","order")), input$Frule)
    })
    output$filter_info <- renderUI({
      fd <- filtered_r(); req(fd)
      ui_filter_info(fd$mv_removed_cols, input$Frule)
    })
    
    params_r <- reactive({
      sel <- selections_r()
      list(
        sample_col = sel$sample, batch_col = sel$batch,
        class_col  = sel$class,  order_col = sel$order,
        withheld_cols = withheld_r(),
        n_withheld = input$n_withheld %||% 0,
        Frule = input$Frule
      )
    })
    
    observeEvent(input$next_correction, {
      validate(need(!is.null(cleaned_r()), "Missing cleaned data"),
               need(!is.null(filtered_r()), "Missing filtered data"))
      updateTabsetPanel(session$rootScope(), "main_steps", "tab_correct")
    })
    
    # module output
    list(
      cleaned  = cleaned_r,
      filtered = filtered_r,
      params   = params_r
    )
  })
}

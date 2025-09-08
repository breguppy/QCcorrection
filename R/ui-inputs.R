#' File upload
#' @keywords internal
#' @noRd
ui_file_upload <- function(ns) {
    fileInput(
      ns("file1"),
      "Choose Raw Data File (.csv, .xls, or .xlsx)",
      accept = c(".csv", ".xls", ".xlsx"),
      buttonLabel = "Browse...",
      placeholder = "No file selected"
    )
}

#' Reusable titled sidebar
#' @keywords internal
#' @noRd
ui_sidebar_block <- function(title, ..., help = NULL, width = 400) {
  bslib::sidebar(
    tags$h4(title),
    ...,
    if (!is.null(help)) lapply(help, tags$h6),
    width = width
  )
}

#' Scrollable table wrapper
#' @keywords internal
#' @noRd
ui_table_scroll <- function(outputId, ns, height = "400px") {
  div(
    style = paste0("overflow:auto; max-height:", height, ";"),
    tableOutput(ns(outputId))
  )
}

#' Column selection for meta data
#' @keywords internal
#' @noRd
ui_nonmet_cols <- function(cols, ns = identity) {
  dropdown_choices <- c("Select a column..." = "", cols)
  
  tagList(
    tooltip(
      selectInput(ns("sample_col"), "sample column", dropdown_choices, ""),
      "Column that contains unique sample names.",
      placement ="right"
    ),
    tooltip(
      selectInput(ns("batch_col"), "batch column", dropdown_choices, ""),
      "Column that contains batch information.",
      placement ="right"
    ),
    tooltip(
      selectInput(ns("class_col"), "class column", dropdown_choices, ""),
      "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.",
      placement ="right"
    ),
    tooltip(
      selectInput(
        ns("order_col"),
        "injection order column",
        dropdown_choices,
        ""
      ),
      "Column that indicates injection order.",
      placement = "right"
    )
  )
}

#' Toggle for withholding extra columns
#' @keywords internal
#' @noRd
ui_withhold_toggle <- function(ns) {
  tooltip(
    checkboxInput(ns("withhold_cols"),
                  "Withhold additional columns from correction?", FALSE),
    "Select if there are extra non-metabolite or specific metabolite columns to withhold.",
    "right"
  )
}

#' Count input for how many columns to withhold
#' @keywords internal
#' @noRd
ui_withhold_count <- function(ns, max_withhold) {
  numericInput(ns("n_withhold"),
               "Number of columns to withhold",
               value = if (max_withhold > 0) 1 else 0,
               min   = 0,
               max   = max_withhold,
               step  = 1)
}

#' missing value filter slider
#' @keywords internal
#' @noRd
ui_filter_slider <- function(ns) {
  tooltip(
    sliderInput(ns("Frule"), "Acceptable % missing per metabolite", 0, 100, 20),
    "Metabolites above this missing % are removed.", "right"
  )
}

#' Repeated selectors for which columns to withhold
#' @param ids character vector of input ids to render (e.g., "withhold_col_1")
#' @param cols candidate column names
#' @param prev named character of previous selections for each id (same length as ids)
#' @keywords internal
#' @noRd
ui_withhold_selectors <- function(ids, cols, prev, ns) {
  if (!length(ids)) return(NULL)
  # keep uniqueness across the repeated selects
  lapply(seq_along(ids), function(i) {
    id    <- ids[i]
    prior <- prev[[i]] %||% ""
    other <- setdiff(prev, prior)
    choices_i <- c("Select a column..." = "", setdiff(cols, other))
    selectInput(
      ns(id),
      label   = paste("Select column to withhold #", i),
      choices = choices_i,
      selected = if (nzchar(prior) && prior %in% choices_i) prior else ""
    )
  })
}

#' Impute missing QC value options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_qc_impute <- function(df, metab_cols, ns = identity) {
  qc_df <- df %>% filter(df$class == "QC")
  has_qc_na <- any(is.na(qc_df[, metab_cols]))
  
  if (has_qc_na) {
    radioButtons(
      ns("qcImputeM"),
      "QC Imputation Method",
      list(
        "metabolite median" = "median",
        "metabolite mean" = "mean",
        "QC-metabolite median" = "class_median",
        "QC-metabolite mean" = "class_mean",
        "minimum value" = "min",
        "half minimum value" = "minHalf",
        "KNN" = "KNN",
        "zero" = "zero"
      ),
      "median",
      FALSE
    )
  } else {
    tags$div(icon("check-circle", class = "text-success"),
             span("No QC missing values"))
  }
}

#' Impute missing sample value options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_sample_impute <- function(df, metab_cols, ns = identity) {
  sam_df <- df %>% filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))
  
  if (has_sam_na) {
    if (num_classes > 1) {
      radioButtons(
        inputId = ns("samImputeM"),
        label = "Sample Imputation Method",
        choices = list(
          "metabolite median" = "median",
          "metabolite mean" = "mean",
          "class-metabolite median" = "class_median",
          "class-metabolite mean" = "class_mean",
          "minimum value" = "min",
          "half minimum value" = "minHalf",
          "KNN" = "KNN",
          "zero" = "zero"
        ),
        selected = "median",
        inline = FALSE
      )
    } else {
      radioButtons(
        inputId = ns("samImputeM"),
        label = "Sample Imputation Method",
        choices = list(
          "metabolite median" = "median",
          "metabolite mean" = "mean",
          "minimum value" = "min",
          "half minimum value" = "minHalf",
          "KNN" = "KNN",
          "zero" = "zero"
        ),
        selected = "median",
        inline = FALSE
      )
    }
  } else {
    tags$div(icon("check-circle", class = "text-success"),
             span("No Sample missing values"))
  }
}

#' Correction method selection options for section 2.1 Choose Correction settings
#' @keywords internal
#' @noRd
ui_correction_method <- function(df, ns = identity) {
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))
  
  if (num_batches == 1) {
    if (any(qc_per_batch$qc_in_batch <= 5)) {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list("Local Polynomial Fit (LOESS)" = "LOESS"),
        selected = "LOESS"
      )
    } else {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS"
        ),
        selected = "RF"
      )
    }
  } else {
    if (any(qc_per_batch$qc_in_batch < 5)) {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS"
        ),
        selected = "RF"
      )
    } else {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Random Forest" = "RF",
          "Local Polynomial Fit (LOESS)" = "LOESS",
          "Batchwise Random Forest" = "BW_RF",
          "Batchwise Local polynomial fit (LOESS)" = "BW_LOESS"
        ),
        selected = "RF"
      )
    }
  }
}
#' @keywords internal
#' @noRd

# Non-metabolite column selection for meta data.
ui_nonmet_cols <- function(cols, ns = identity) {
  dropdown_choices <- c("Select a column..." = "", cols)
  
  tagList(
    tooltip(
      selectInput(ns("sample_col"), "sample column", dropdown_choices, ""),
      "Column that contains unique sample names.",
      "right"
    ),
    tooltip(
      selectInput(ns("batch_col"), "batch column", dropdown_choices, ""),
      "Column that contains batch information.",
      "right"
    ),
    tooltip(
      selectInput(ns("class_col"), "class column", dropdown_choices, ""),
      "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.",
      "right"
    ),
    tooltip(
      selectInput(
        ns("order_col"),
        "injection order column",
        dropdown_choices,
        ""
      ),
      "Column that indicates injection order.",
      "right"
    )
  )
}

# Impute missing QC value options for section 2.1 Choose Correction settings
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

# Impute missing sample value options for section 2.1 Choose Correction settings
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

# Correction method selection options for section 2.1 Choose Correction settings
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
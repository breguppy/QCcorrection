# Helper function for app.R
library(shiny)
library(bslib)
library(dplyr)
library(shinycssloaders)
library(purrr)

#–– UI snippets ––#
# Non-metabolite column selection for meta data.
nonMetColSelectionUI <- function(cols, ns = identity) {
  dropdown_choices <- c("Select a column..." = "", cols)
  
  tagList(
    tooltip(
      selectInput(ns("sample_col"), "sample column", dropdown_choices, ""),
      "Column that contains unique sample names.", "right"
    ),
    tooltip(
      selectInput(ns("batch_col"), "batch column", dropdown_choices, ""),
      "Column that contains batch information.", "right"
    ),
    tooltip(
      selectInput(ns("class_col"), "class column", dropdown_choices, ""),
      "Column that indicates the type of sample. Must contain QC samples labeled as 'NA', 'QC', 'Qc', or 'qc'.", "right"
    ),
    tooltip(
      selectInput(ns("order_col"), "injection order column", dropdown_choices, ""),
      "Column that indicates injection order.", "right"
    )
  )
}

# Column-warning text for section 1.2 Select non-metabolite columns
columnWarningUI <- function(data, selected) {
  warnings <- list()
  
  if (any(selected == "")) {
    warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkorange; font-weight:bold;",
                                                  icon("exclamation-triangle"),
                                                  " Please select all four columns.")
  } else if (length(unique(selected)) < 4) {
    warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold;",
                                                  icon("exclamation-triangle"),
                                                  " Each selected column must be unique.")
  }
  
  if (!any(selected == "") && length(unique(selected)) == 4) {
    samp_vec  <- data[[selected[1]]]
    order_vec <- data[[selected[4]]]
    
    if (anyDuplicated(samp_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
                                                    icon("exclamation-triangle"),
                                                    " Duplicate sample names detected!")
    }
    
    if (anyDuplicated(order_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
                                                    icon("exclamation-triangle"),
                                                    " Duplicate order values detected!")
    }
  }
  
  if (length(warnings) == 0) {
    return(NULL)
  } else if (length(warnings) == 1) {
    return(warnings[[1]])
  } else {
    return(do.call(tagList, warnings))
  }
}

# Metric card for 1.2 Select non-metabolite columns
metric_card <- function(label, value) {
  div(
    style = "background:#f8f9fa; padding:10px; border-radius:8px; flex:1;",
    p(style = "font-size:1.4em; font-weight:bold; margin:0;", value),
    h5(style = "margin:0;", label)
  )
}

# Basic info for section 1.2 Select non-metabolite columns
basicInfoUI <- function(df, replacement_counts) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  n_metab = length(metab_cols)
  n_missv = sum(is.na(df[, metab_cols]))
  n_qcs   = sum(df$class == "QC")
  n_samp  = sum(df$class != "QC")
  n_bat   = n_distinct(df$batch)
  n_class = n_distinct(df$class[df$class != "QC"])
  class_list <- sort(unique(df$class[df$class != "QC"]))
  perc_missv <- round(100 * (n_missv / ((n_samp + n_qcs) * n_metab)), digits = 2)
  
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_class = sum(class == "QC"), .groups = "drop")
  
  total_replaced <- sum(replacement_counts$non_numeric_replaced +
                          replacement_counts$zero_replaced)
  
  class_badges <- tags$div(style = "display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;", lapply(class_list, function(cls) {
    tags$span(style = "background-color: #e9ecef; padding: 5px 10px; border-radius: 12px;", as.character(cls))
  }))
  
  na_table_html <- tags$div(
    style = "margin-top: 15px;",
    tags$h5("Number of QC Samples per Batch"),
    tags$table(class = "table table-bordered table-sm", tags$thead(tags$tr(
      tags$th("Batch"), tags$th("QCs in Batch")
    )), tags$tbody(lapply(1:nrow(qc_per_batch), function(i) {
      tags$tr(tags$td(as.character(qc_per_batch$batch[i])),
              tags$td(qc_per_batch$qc_in_class[i]))
    })))
  )
  
  tagList(
    if (total_replaced > 0) {
      tags$span(
        style = "color: darkred; font-weight: bold;",
        paste(
          total_replaced,
          "non-numeric or zero metabolite values are counted as missing."
        )
      )
    },
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 10px;",
      # left side
      tags$div(
        style = "display: grid; grid-template-columns: repeat(1, 1fr); gap: 20px;",
        
        # top left
        tags$div(
          style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-top 15px;",
          metric_card("Metabolite Columns", n_metab),
          metric_card("Missing Values", paste0(n_missv, " (", perc_missv, "%)")),
          metric_card("QC Samples", n_qcs),
          metric_card("Samples", n_samp),
          metric_card("Batches", n_bat),
          metric_card("Classes", n_class),
        ),
        
        #bottom left
        tags$div(style = "flex: 1; min-width: 250px;", tags$h5("Unique Classes"), class_badges)
        
      ),
      
      # right half
      # QC per batch table column
      tags$div(
        style = "flex: 1; min-width: 250px;",
        tags$h5("Number of QC Samples per Batch"),
        tags$table(class = "table table-bordered table-sm", tags$thead(tags$tr(
          tags$th("Batch"), tags$th("QCs in Batch")
        )), tags$tbody(lapply(1:nrow(qc_per_batch), function(i) {
          tags$tr(tags$td(as.character(qc_per_batch$batch[i])),
                  tags$td(qc_per_batch$qc_in_class[i]))
        })))
      )
    )
  )
}

# Filter info for section 1.4 Filter Raw Data
filterInfoUI <- function(mv_removed, Frule) {
  if (length(mv_removed) == 0) {
    tags$div(
      style = "flex: 1; padding-right: 10px;",
      tags$span(
        style = "color:darkgreen;font-weight:bold;",
        paste0("No metabolites removed for missing value percentage above ",
              Frule, "%.")
      )
    )
  } else {
    tags$div(
      style = "flex: 1; padding-right: 10px;",
      tags$span(
        style = "color:darkorange;font-weight:bold;",
        paste0(
          length(mv_removed),
          " metabolites removed based on missing value percentage above ",
          Frule, "%"
        )
      ),
      tags$ul(lapply(mv_removed, tags$li))
    )
  }
  
}

# QC missing value warning for section 2.1 Choose Correction settings
qcMissingValueWarning <- function(df) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  qc_idx <- which(df$class == "QC")
  n_missv = sum(is.na(df[qc_idx, metab_cols]))
  
  if (n_missv > 0) {
    tags$span(
      style = "color:darkred; font-weight:bold;",
      icon("exclamation-triangle"),
      paste(" ", n_missv, " values missing from QC samples")
    )
  } else {
    NULL
  }
  
}

# Impute missing QC value options for section 2.1 Choose Correction settings
qcImputeUI <- function(df, metab_cols, ns = identity) {
  qc_df <- df %>% filter(df$class == "QC")
  has_qc_na <- any(is.na(qc_df[, metab_cols]))
  
  if (has_qc_na) {
    radioButtons(ns("qcImputeM"), "QC Imputation Method", list(
        "metabolite median" = "median",
        "metabolite mean" = "mean",
        "QC-metabolite median" = "class_median",
        "QC-metabolite mean" = "class_mean",
        "minimum value" = "min",
        "half minimum value" = "minHalf",
        "KNN" = "KNN",
        "zero" = "zero"
      ),
      "median", FALSE
    )
  } else {
    tags$div(icon("check-circle", class = "text-success"),
             span("No QC missing values"))
  }
}

# Impute missing sample value options for section 2.1 Choose Correction settings
sampleImputeUI <- function(df, metab_cols, ns =identity) {
  sam_df <- df %>% filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))
  
  if (has_sam_na) {
    if (num_classes > 1){
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
correctionMethodUI <- function(df, ns = identity) {
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))
  
  if (num_batches == 1) {
    if (any(qc_per_batch$qc_in_batch <= 5)) {
      radioButtons(
        inputId = ns("corMethod"),
        label = "Method",
        choices = list(
          "Local Polynomial Fit (LOESS)" = "LOESS"
        ),
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

# Unavailable correction option description for section 2.1 Choose Correction settings
unavailableOptionsUI <- function(df, metab_cols) {
  # If there is only 1 class: class/metabolite impute for samples not available
  # If there is only 1 batch: no batchwise options
  # If there is less than 5 QCs per batch: no batchwise options
  # If there is only 1 batch and less than 5 QCs no RF option
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_batch = sum(class == "QC"), .groups = "drop")
  num_batches <- length(unique(df$batch))
  
  sam_df <- df %>% filter(df$class != "QC")
  has_sam_na <- any(is.na(sam_df[, metab_cols]))
  num_classes <- length(unique(sam_df$class))
  
  unavail_opts <- list()
  if (has_sam_na & (num_classes == 1)) {
    unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Sample Imputation Methods:")
    unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
      icon("circle-xmark", class = "text-danger-emphasis"),
      " class-metabolite median requires more than 1 class."
    )
    unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
    unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
      icon("circle-xmark", class = "text-danger-emphasis"),
      " class-metabolite mean requires more than 1 class."
    )
  }
  if (num_batches == 1) {
    if (any(qc_per_batch$qc_in_batch <= 5)) {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Random Forest requires more QC samples."
      )
      unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise options (Random Forest and LOESS) require more than 1 batch."
      )
    } else {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise options (Random Forest and LOESS) require more than 1 batch."
      )
    }
  } else {
    if (any(qc_per_batch$qc_in_batch < 5)) {
      unavail_opts[[length(unavail_opts) + 1]] <- tags$h6("Unavailable Correction Methods:")
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise Random Forest requires at least 5 QCs per batch."
      )
      unavail_opts[[length(unavail_opts) + 1]] <- tags$br()
      unavail_opts[[length(unavail_opts) + 1]] <- tags$span(
        icon("circle-xmark", class = "text-danger-emphasis"),
        " Batchwise LOESS requires at least 5 QCs per batch."
      )
    }
  }
  
  if (length(unavail_opts) == 0) {
    return(tags$span("All methods available"))
  } else if (length(unavail_opts) == 1) {
    return(unavail_opts[[1]])
  } else {
    return(do.call(tagList, unavail_opts))
  }
}

# Post-correction filtering info for section 2.2 Post-Correction Filtering
postCorFilterInfoUI <- function(filtered_corrected_result, rsd_filter, post_cor_filter) {
  n_removed <- length(filtered_corrected_result$removed_metabolites)
  if (post_cor_filter == FALSE) {
    ui <- list(
      tags$span(
        style = "color: darkorange; font-weight: bold;",
        paste0(
          n_removed,
          " metabolites removed based on QC RSD above ",
          rsd_filter, "%"
        )
      ),
      tags$br(),
      tags$ul(
        lapply(filtered_corrected_result$removed_metabolites, function(name) {
          tags$li(name)
        })
      )
    )
  } else {
    ui <- list(
      tags$span(
        style = "color: darkgreen; font-weight: bold;",
        "Metabolites are not filtered by QC RSD."
      )
    )
  }
  do.call(tagList, ui)
}

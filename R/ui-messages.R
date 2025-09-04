#' @keywords internal
#' @noRd

# Column-warning text for section 1.2 Select non-metabolite columns
ui_column_warning <- function(data, selected) {
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

# QC missing value warning for section 2.1 Choose Correction settings
ui_qc_missing_warning <- function(df) {
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

# Unavailable correction option description for section 2.1 Choose Correction settings
ui_unavailable_options <- function(df, metab_cols) {
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
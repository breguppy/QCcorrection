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

ui_rsd_stats <- function(p, d) {
  df_before <- d$filtered$df
  # Determine df_after based on rsd_compare selected by user.
  if (p$rsd_compare == "filtered_cor_data") {
    df_after <- d$filtered_corrected$df
  } else {
    df_after <- d$transformed$df
  }
  if (p$rsd_cal == "met") {
    rsdBefore <- metabolite_rsd(df_before)
    rsdAfter <- metabolite_rsd(df_after)
  } else {
    rsdBefore <- class_metabolite_rsd(df_before)
    rsdAfter <- class_metabolite_rsd(df_after)
  }
  rsd_stats <- delta_rsd_stats(rsdBefore, rsdAfter)
  
  df <- data.frame(
    Metric = c("Average &Delta; QC RSD", 
               "Median &Delta; QC RSD", 
               "Average &Delta; Sample RSD", 
               "Median &Delta; Sample RSD"),
    Value  = c(rsd_stats$avg_delta_qc,
               rsd_stats$med_delta_qc,
               rsd_stats$avg_delta_sample,
               rsd_stats$med_delta_sample)
  )
  
  df$Value <- sprintf("%.2f%%", df$Value)
  
  htmltools::tagList(
    htmltools::tags$table(
      style = "border-collapse: collapse; margin-top:10px;",
      htmltools::tags$thead(
        htmltools::tags$tr(
          htmltools::tags$th("Metric",  style="padding:4px 12px; text-align:left; border-bottom:1px solid #ccc;"),
          htmltools::tags$th("Value",   style="padding:4px 12px; text-align:right; border-bottom:1px solid #ccc;")
        )
      ),
      htmltools::tags$tbody(
        lapply(seq_len(nrow(df)), function(i) {
          htmltools::tags$tr(
            htmltools::tags$td(HTML(df$Metric[i]), style="padding:4px 12px; text-align:left;"),
            htmltools::tags$td(df$Value[i],  style="padding:4px 12px; text-align:right;")
          )
        })
      )
    )
  )
}

ui_md_outlier <- function(p, d, digits = 2, max_contrib = 5) {
  if (p$pca_compare == "filtered_cor_data") {
    df <- d$filtered_corrected$df
    after <- d$filtered_corrected
    compared_to <- "Correction"
  } else {
    df <- d$transformed$df
    after <- d$transformed
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(df), c("sample", "batch", "class", "order"))
  shiny::validate(
    shiny::need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    shiny::need(nrow(df) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- df[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  shiny::validate(need(
    any(keep),
    "All metabolite columns are constant/invalid after filtering."
  ))
  
  # before data cannot have any missing values.
  #before <- d$imputed
  md_res <- md_outliers_by_group(p, after, after, robust = TRUE, alpha = 0.95)
  stopifnot(all(c("dataset","group","sample","MD2","is_outlier","top_contributors") %in% names(md_res)))
  
  df <- md_res |>
    dplyr::filter(isTRUE(.data$is_outlier), is.finite(.data$MD2)) |>
    dplyr::mutate(
      Sample       = .data$sample,
      Mahalanobis  = round(sqrt(.data$MD2), digits),
      Contributors = vapply(.data$top_contributors,
                            function(x) paste(utils::head(x %||% character(), max_contrib), collapse = ", "),
                            character(1))
    ) |>
    dplyr::select(.data$dataset, .data$group, .data$Sample, .data$Mahalanobis, .data$Contributors)
  
  if (nrow(df) == 0L) {
    return(htmltools::tags$div("No outlier samples detected in PCA."))
  }
  
  by_block <- split(df, list(df$dataset, df$group), drop = TRUE)
  
  make_table <- function(block) {
    ds <- unique(block$dataset)
    gp <- unique(block$group)
    htmltools::tags$div(
      htmltools::tags$h4(sprintf("%s — %s", ds, gp)),
      htmltools::tags$table(
        class = "table table-sm table-striped",
        style = "width:100%; table-layout:fixed; word-wrap:break-word;",
        htmltools::tags$thead(
          htmltools::tags$tr(
            htmltools::tags$th("Sample"),
            htmltools::tags$th("Mahalanobis distance"),
            htmltools::tags$th("Top contributors")
          )
        ),
        htmltools::tags$tbody(
          lapply(seq_len(nrow(block)), function(i) {
            htmltools::tags$tr(
              htmltools::tags$td(htmltools::htmlEscape(block$Sample[i])),
              htmltools::tags$td(format(block$Mahalanobis[i], nsmall = digits, trim = TRUE)),
              htmltools::tags$td(htmltools::htmlEscape(block$Contributors[i]))
            )
          })
        )
      )
    )
  }
  
  htmltools::tagList(lapply(by_block, make_table))
}
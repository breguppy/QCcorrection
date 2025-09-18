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

ui_outliers <- function(p, d, confirmations = NULL, sample_md = NULL,
                        max_top = 5, digits = 2) {
  if (p$out_data == "filtered_cor_data") {
    df <- d$filtered_corrected$df
  } else {
    df <- d$transformed$df
  }
  res <- detect_qc_aware_outliers(df, group_nonqc_by_class = p$sample_grouping)
  
  if (!requireNamespace("htmltools", quietly = TRUE)) stop("htmltools is required")
  tags <- htmltools::tags
  
  # Accept either the full result list, or separate pieces
  if (is.null(res)) res <- list(confirmations = confirmations, sample_md = sample_md)
  if (!is.list(res)) stop("ui_outliers: supply the full result list or confirmations+sample_md")
  
  conf <- res$confirmations
  md   <- res$sample_md
  if (is.null(conf) || is.null(md)) stop("ui_outliers: missing confirmations or sample_md")
  
  # choose grouping column present in md (group_id for new code, else class)
  group_col <- if ("group_id" %in% names(md)) "group_id" else "class"
  
  conf_ok <- subset(conf, decision == "confirm")
  if (nrow(conf_ok) == 0) return(tags$span("No outliers detected"))
  
  conf_ok$abs_z <- abs(conf_ok$z)
  
  top_by_sample <- lapply(split(conf_ok, conf_ok$sample), function(d) {
    d <- d[order(-d$abs_z), ]
    d <- head(d, max_top)
    data.frame(
      sample   = d$sample[1],
      groupval = d[[group_col]][1],
      top_mets = paste0(
        d$metabolite,
        " (z=", sprintf("%.*f", digits, d$z),
        ", QC RSD=", ifelse(is.na(d$qc_rsd), "NA", sprintf("%.*f%%", digits, d$qc_rsd)),
        ", ", d$method, ":p=", ifelse(is.na(d$p_value), "NA", sprintf("%.*f", digits, d$p_value)),
        ")",
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  })
  top_by_sample <- do.call(rbind, top_by_sample)
  
  md_keep <- md[, c("sample", group_col, "md", "cutoff", "flagged")]
  names(md_keep)[names(md_keep) == group_col] <- "groupval"
  out_tbl <- merge(top_by_sample, md_keep, by = c("sample","groupval"), all.x = TRUE)
  
  out_tbl$flagged <- with(out_tbl, isTRUE(flagged) | (!is.na(md) & md > cutoff))
  out_tbl <- out_tbl[order(-as.numeric(out_tbl$flagged), -out_tbl$md), ]
  
  header <- tags$thead(
    tags$tr(
      tags$th("Sample"), tags$th("Group"),
      tags$th("Mahalanobis"), tags$th("Cutoff"),
      tags$th("Flagged"), tags$th(paste0("Top metabolites (max ", max_top, ")"))
    )
  )
  body_rows <- apply(out_tbl, 1, function(r) {
    tags$tr(
      tags$td(r[["sample"]]),
      tags$td(r[["groupval"]]),
      tags$td(ifelse(is.na(r[["md"]]), "NA", sprintf("%.*f", digits, as.numeric(r[["md"]])))),
      tags$td(ifelse(is.na(r[["cutoff"]]), "NA", sprintf("%.*f", digits, as.numeric(r[["cutoff"]])))),
      tags$td(ifelse(isTRUE(as.logical(r[["flagged"]])), "yes", "no")),
      tags$td(r[["top_mets"]])
    )
  })
  
  tags$table(
    style = "border-collapse:collapse; width:100%; font-size:90%;",
    tags$style(htmltools::HTML("
      table, th, td { border: 1px solid #ccc; }
      th, td { padding: 6px 8px; vertical-align: top; }
      th { background:#f7f7f7; text-align:left; }
    ")),
    header, tags$tbody(body_rows)
  )
}

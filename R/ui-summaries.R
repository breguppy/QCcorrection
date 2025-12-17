#' @keywords internal
#' @noRd

# Metric card for 1.2 Select non-metabolite columns
metric_card <- function(label, value) {
  div(
    style = "background:#f8f9fa; padding:10px; border-radius:8px; flex:1;",
    p(style = "font-size:1.4em; font-weight:bold; margin:0;", value),
    h5(style = "margin:0;", label)
  )
}
# Reusable warning card generator
warn_card <- function(title, body, body_tags = NULL) {
  tags$div(
    class = "card border-warning mb-3",
    style = "margin-top: 10px;",
    tags$div(class = "card-header", title),
    tags$div(
      class = "card-body",
      tags$p(class = "card-text", body),
      body_tags
    )
  )
}

# Basic info for section 1.2 Select non-metabolite columns
ui_basic_info <- function(df,
                          replacement_counts,
                          non_numeric_cols,
                          duplicate_mets = NULL,
                          high_corr_mets = NULL) {
  
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  n_metab    <- length(metab_cols)
  n_missv    <- sum(is.na(df[, metab_cols]))
  n_qcs      <- sum(df$class == "QC")
  n_samp     <- sum(df$class != "QC")
  n_bat      <- dplyr::n_distinct(df$batch)
  n_class    <- dplyr::n_distinct(df$class[df$class != "QC"])
  class_list <- sort(unique(df$class[df$class != "QC"]))
  perc_missv <- round(100 * (n_missv / ((n_samp + n_qcs) * n_metab)), digits = 2)
  
  qc_per_batch <- df %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(qc_in_class = sum(class == "QC"), .groups = "drop")
  
  total_replaced <- sum(
    replacement_counts$non_numeric_replaced +
      replacement_counts$zero_replaced
  )
  
  class_badges <- tags$div(
    style = "display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;",
    lapply(class_list, function(cls) {
      tags$span(
        style = "background-color: #e9ecef; padding: 5px 10px; border-radius: 12px;",
        as.character(cls)
      )
    })
  )
  
  # ---------- Warning box 1: replaced values ----------
  replaced_card <- NULL
  if (total_replaced > 0) {
    replaced_card <- warn_card(
      title = "Replaced non-numeric or zero metabolite values",
      body  = paste0(
        total_replaced,
        " values were converted to missing (NA) prior to processing."
      )
    )
  }
  
  # ---------- Warning box 2: non-numeric columns ----------
  nonnum_card <- NULL
  if (length(non_numeric_cols) > 0) {
    nonnum_card <- warn_card(
      title = "Non-numerical columns detected",
      body  = "These columns contain all non-numeric values and will be removed prior to processing.", #If you do not want these columns removed, check the box 'withhold additional columns from correction' on the left and add them to list below the box.",
      body_tags = tags$p(
        style = "font-weight: 600; margin-top: 8px;",
        paste(sort(unique(non_numeric_cols)), collapse = ", ")
      )
    )
  }
  
  # ---------- Warning box 3: duplicate metabolites ----------
  duplicate_card <- NULL
  if (!is.null(duplicate_mets) && nrow(duplicate_mets) > 0) {
    
    dup_badges <- tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px;",
      lapply(seq_len(nrow(duplicate_mets)), function(i) {
        pair <- duplicate_mets[i, , drop = FALSE]
        tags$span(
          style = paste(
            "background-color: #ff8989;",
            "border: 1px solid #ff8989;",
            "padding: 4px 8px;",
            "border-radius: 12px;",
            "font-size: 0.85rem;"
          ),
          sprintf("%s \u2248 %s", pair$col1, pair$col2)
        )
      })
    )
    
    duplicate_card <- warn_card(
      title = "Potential duplicate metabolites",
      body  = sprintf(
        "%d column pairs appear equal or nearly equal based on non-missing values.",
        nrow(duplicate_mets)
      ),
      body_tags = dup_badges
    )
  }
  
  
  # ---------------------------------------------------
  # MAIN UI SECTION
  # ---------------------------------------------------
  tagList(
    replaced_card,
    nonnum_card,
    duplicate_card,
    
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 10px;",
      
      # left section
      tags$div(
        style = "display: grid; grid-template-columns: repeat(1, 1fr); gap: 20px;",
        
        # metrics grid
        tags$div(
          style = "display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-top: 15px;",
          metric_card("Metabolite Columns", n_metab),
          metric_card("Missing Values", paste0(n_missv, " (", perc_missv, "%)")),
          metric_card("QC Samples", n_qcs),
          metric_card("Samples", n_samp),
          metric_card("Batches", n_bat),
          metric_card("Classes", n_class)
        ),
        
        # class list badges
        tags$div(
          style = "flex: 1; min-width: 250px;",
          tags$h5("Unique Classes"),
          class_badges
        )
      ),
      
      # right section
      tags$div(
        style = "flex: 1; min-width: 250px;",
        tags$h5("Number of QC Samples per Batch"),
        tags$table(
          class = "table table-bordered table-sm",
          tags$thead(
            tags$tr(tags$th("Batch"), tags$th("QCs in Batch"))
          ),
          tags$tbody(
            lapply(seq_len(nrow(qc_per_batch)), function(i) {
              tags$tr(
                tags$td(as.character(qc_per_batch$batch[i])),
                tags$td(qc_per_batch$qc_in_class[i])
              )
            })
          )
        )
      )
    )
  )
}



# Filter info for section 1.4 Filter Raw Data
ui_filter_info <- function(mv_removed, mv_cutoff, qc_missing_mets) {
  left_col <- if (length(mv_removed) == 0) {
    tags$div(style = "flex: 1; padding-right: 10px;",
             tags$span(
               style = "color:darkgreen;font-weight:bold;",
               paste0(
                 "No metabolites removed for missing value percentage above ",
                 mv_cutoff,
                 "%."
               )
             ))
  } else {
    tags$div(
      style = "flex: 1; padding-right: 10px;",
      tags$span(
        style = "color:darkorange;font-weight:bold;",
        paste0(
          length(mv_removed),
          " metabolite(s) removed based on missing value percentage above ",
          mv_cutoff,
          "%"
        )
      ),
      tags$ul(lapply(mv_removed, tags$li))
    )
  }
  
  right_col <- if(length(qc_missing_mets) == 0) {
    tags$div(
      style = "flex:1; padding-left:10px;",
      tags$span(style = "color:darkgreen; font-weight:bold;",
                "No metabolites have missing values in QC samples after filtering."))
  } else {
    tags$div(
      style = "flex:1; padding-left:10px;",
      tags$span(style = "color:darkorange; font-weight:bold;",
                paste0(length(qc_missing_mets),
                       " metabolite(s) with at least one QC missing value after filtering.")),
      tags$ul(lapply(qc_missing_mets, tags$li)))
  }
  
  tags$div(
    style = "display:flex; gap:16px; align-items:flex-start;",
    left_col, right_col
  )
}


ui_corr_range_info <- function(all_cor, range) {
  high_corr_mets <- filter_correlation_pairs_by_range(all_cor, range)
  correlated_card <- tags$div(style = "flex: 1; padding-right: 10px;",
                              tags$span(
                                style = "color:darkgreen;font-weight:bold;",
                                "No metabolite pairs within this correlation range."
                              ))
  if (!is.null(high_corr_mets) && nrow(high_corr_mets) > 0) {
    
    cor_badges <- tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px;",
      lapply(seq_len(nrow(high_corr_mets)), function(i) {
        pair <- high_corr_mets[i, , drop = FALSE]
        tags$span(
          style = paste(
            "background-color: #fff3cd;",
            "border: 1px solid #ffeeba;",
            "padding: 4px 8px;",
            "border-radius: 12px;",
            "font-size: 0.85rem;"
          ),
          sprintf("%s \u221D %s", pair$col1, pair$col2)
        )
      })
    )
    
    correlated_card <- warn_card(
      title = "Correlated Metabolites",
      body  = sprintf(
        "%d column pairs have correlation (Pearson's r) within the selected range.",
        nrow(high_corr_mets)
      ),
      body_tags = cor_badges
    )
  
  } 
  correlated_card
}

# Post-correction filtering info for section 2.2 Post-Correction Filtering
ui_postcor_filter_info <- function(filtered_corrected_result,
                                   remove_imputed,
                                   rsd_filter,
                                   post_cor_filter) {
  if (isTRUE(remove_imputed)) {
    removed <- filtered_corrected_result$removed_metabolites_mv
  } else {
    removed <- filtered_corrected_result$removed_metabolites_no_mv
  }
  n_removed <- length(removed)
  
  # get ISTD/ITSD metabolites
  is_istd <- grepl("^(ISTD|ITSD)", removed, ignore.case = TRUE)
  istd_names <- removed[is_istd]
  n_istd <- length(istd_names)
  
  # optional warning banner
  warning_ui <- NULL
  if (n_istd > 0) {
    warning_ui <- tags$div(
      class = "alert alert-danger",
      style = "margin-bottom: 10px;",
      tags$strong(paste0(n_istd, " internal standard(s) with QC RSD above ", rsd_filter, "%: ")),
      tags$ul(
        lapply(istd_names, tags$li)
      )
    )
  }
  
  if (post_cor_filter == FALSE) {
    ui <- list(
      warning_ui,
      tags$span(
        style = "color: darkorange; font-weight: bold;",
        paste0(
          n_removed,
          " metabolite(s) removed based on QC RSD above ",
          rsd_filter,
          "%"
        )
      ),
      tags$br(),
      tags$ul(
        lapply(removed, tags$li)
      )
    )
  } else {
    ui <- list(
      warning_ui,
      tags$span(
        style = "color: darkgreen; font-weight: bold;",
        "Metabolites are not filtered by QC RSD."
      )
    )
  }
  
  do.call(tagList, ui)
}

ui_outliers <- function(p, d, 
                        top_n        = 10L,
                        sample_col   = "sample",
                        class_col    = "class",
                        digits_z     = 2L,
                        digits_T2    = 2L,
                        pca_output_id = "hotelling_pca",
                        ns           = identity) {
  df <- if (p$out_data == "filtered_cor_data") {
      d$filtered_corrected$df_no_mv
  } else {
    d$transformed$df_no_mv
  }
  
  detect_result <- detect_hotelling_nonqc_dual_z(df, p)
  
  # Extract extreme values and full data
  ev   <- detect_result$extreme_values
  dres <- detect_result$data
  
  if (is.null(ev)) {
    stop("detect_result$extreme_values is NULL. Did you pass the correct object?")
  }
  
  # Metric values
  n_outlier_samples  <- sum(dres$is_outlier_sample, na.rm = TRUE)
  n_extreme_values   <- nrow(ev)
  
  # Metric cards
  cards <- shiny::div(
    style = "display:flex; gap:10px; margin-bottom:10px;",
    metric_card("Samples outside the Hotelling's T^2 95% limit", n_outlier_samples),
    metric_card("Potential extreme metabolite values", n_extreme_values)
  )
  
  # If no extreme values, just show cards + PCA + message
  if (nrow(ev) == 0L) {
    return(shiny::tagList(
      shiny::plotOutput(ns(pca_output_id), height = "350px"),
      cards,
      shiny::tags$em("No extreme metabolite values detected in outlier samples.")
    ))
  }
  
  # Check required columns from detect_hotelling_nonqc_dual_z()
  required_cols <- c(
    sample_col,
    class_col,
    "metabolite",
    "z_global",
    "abs_z_global",
    "z_class",
    "abs_z_class",
    "T2"
  )
  missing_cols <- setdiff(required_cols, names(ev))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in extreme_values: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Sort by |z_global|, then |z_class|, then T2, descending
  ev_sorted <- ev[order(-ev$abs_z_global, -ev$abs_z_class, -ev$T2), , drop = FALSE]
  
  # Take top_n rows
  ev_top <- head(ev_sorted, top_n)
  
  # Format numeric values
  z_g_fmt    <- formatC(ev_top$z_global,      format = "f", digits = digits_z)
  absz_g_fmt <- formatC(ev_top$abs_z_global,  format = "f", digits = digits_z)
  z_c_fmt    <- formatC(ev_top$z_class,       format = "f", digits = digits_z)
  absz_c_fmt <- formatC(ev_top$abs_z_class,   format = "f", digits = digits_z)
  T2_fmt     <- formatC(ev_top$T2,            format = "f", digits = digits_T2)
  
  # Build table rows
  rows <- lapply(seq_len(nrow(ev_top)), function(i) {
    shiny::tags$tr(
      shiny::tags$td(ev_top[[sample_col]][i]),
      shiny::tags$td(ev_top[[class_col]][i]),
      shiny::tags$td(ev_top$metabolite[i]),
      shiny::tags$td(z_g_fmt[i]),
      shiny::tags$td(absz_g_fmt[i]),
      shiny::tags$td(z_c_fmt[i]),
      shiny::tags$td(absz_c_fmt[i]),
      shiny::tags$td(T2_fmt[i])
    )
  })
  
  # Build full table with header
  table_tag <- shiny::tags$table(
    class = "table table-striped table-condensed table-hover",
    shiny::tags$thead(
      shiny::tags$tr(
        shiny::tags$th("Sample"),
        shiny::tags$th("Class"),
        shiny::tags$th("Metabolite"),
        shiny::tags$th("Global z-score"),
        shiny::tags$th("|z| (global)"),
        shiny::tags$th("Class z-score"),
        shiny::tags$th("|z| (class)"),
        shiny::tags$th("Mahalanobis^2")
      )
    ),
    shiny::tags$tbody(rows)
  )
  
  shiny::tagList(
    # PCA plot placeholder – actual plot comes from renderPlot()
    shiny::plotOutput(ns(pca_output_id), height = "350px"),
    cards,
    shiny::tags$span(
      "Top 10 potential extreme values are listed below. ",
      "The full list of potential extreme values 'extreme_values_*today's_date*.xlsx' ",
      "is available for download."
    ),
    table_tag
  )
}




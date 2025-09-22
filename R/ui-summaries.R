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

# Basic info for section 1.2 Select non-metabolite columns
ui_basic_info <- function(df, replacement_counts) {
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

# Post-correction filtering info for section 2.2 Post-Correction Filtering
ui_postcor_filter_info <- function(filtered_corrected_result,
                                   rsd_filter,
                                   post_cor_filter) {
  n_removed <- length(filtered_corrected_result$removed_metabolites)
  if (post_cor_filter == FALSE) {
    ui <- list(
      tags$span(
        style = "color: darkorange; font-weight: bold;",
        paste0(
          n_removed,
          " metabolites removed based on QC RSD above ",
          rsd_filter,
          "%"
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
      tags$span(style = "color: darkgreen; font-weight: bold;", "Metabolites are not filtered by QC RSD.")
    )
  }
  do.call(tagList, ui)
}

ui_outliers <- function(p, d, confirmations = NULL, sample_md = NULL,
                        digits = 2, n_top = 10) {
  df <- if (p$out_data == "filtered_cor_data") d$filtered_corrected$df else d$transformed$df
  res <- detect_qc_aware_outliers(df, group_nonqc_by_class = p$sample_grouping)
  
  if (!requireNamespace("htmltools", quietly = TRUE)) stop("htmltools is required")
  tags <- htmltools::tags
  
  if (is.null(res)) res <- list(confirmations = confirmations, sample_md = sample_md)
  if (!is.list(res)) stop("ui_outliers: supply the full result list or confirmations+sample_md")
  
  conf <- res$confirmations
  md   <- res$sample_md
  if (is.null(conf) || is.null(md)) stop("ui_outliers: missing confirmations or sample_md")
  
  group_col <- if ("group_id" %in% names(md)) "group_id" else "class"
  conf_ok <- subset(conf, decision == "confirm")
  
  n_samples_with_outlier <- if (nrow(conf_ok)) length(unique(conf_ok$sample)) else 0L
  n_candidate_values     <- nrow(conf_ok)
  
  if (nrow(conf_ok) == 0) {
    return(tagList(
      tags$div(
        style = "display:flex; gap:16px; margin-bottom:12px;",
        metric_card("Samples with at least 1 potential outlier", n_samples_with_outlier),
        metric_card("Candidate outlier values", n_candidate_values)
      ),
      tags$span("No outliers detected")
    ))
  }
  
  # join Mahalanobis info for each (sample, group)
  md_keep <- md[, c("sample", group_col, "md", "cutoff", "flagged")]
  names(md_keep)[names(md_keep) == group_col] <- "groupval"
  
  conf_ok$abs_z <- abs(conf_ok$z)
  conf_ok$groupval <- conf_ok[[group_col]]
  top <- merge(conf_ok, md_keep, by = c("sample","groupval"), all.x = TRUE)
  
  # rank by |z|, then MD desc
  top <- top[order(-top$abs_z, -as.numeric(top$md)), ]
  top <- head(top, n_top)
  
  # format rows
  header <- tags$thead(
    tags$tr(
      tags$th("#"),
      tags$th("Sample"),
      tags$th("Group"),
      tags$th("Metabolite"),
      tags$th("|z|"),
      tags$th("z"),
      tags$th("QC RSD"),
      tags$th("Method"),
      tags$th("p"),
      tags$th("Mahalanobis"),
      tags$th("Cutoff")
    )
  )
  
  body_rows <- lapply(seq_len(nrow(top)), function(i) {
    r <- top[i, ]
    tags$tr(
      tags$td(i),
      tags$td(r$sample),
      tags$td(r$groupval),
      tags$td(r$metabolite),
      tags$td(sprintf("%.*f", digits, r$abs_z)),
      tags$td(sprintf("%.*f", digits, r$z)),
      tags$td(ifelse(is.na(r$qc_rsd), "NA", sprintf("%.*f%%", digits, r$qc_rsd))),
      tags$td(r$method),
      tags$td(ifelse(is.na(r$p_value), "NA", sprintf("%.*f", digits, r$p_value))),
      tags$td(ifelse(is.na(r$md), "NA", sprintf("%.*f", digits, r$md))),
      tags$td(ifelse(is.na(r$cutoff), "NA", sprintf("%.*f", digits, r$cutoff)))
    )
  })
  
  tagList(
    tags$div(
      style = "display:flex; gap:16px; margin-bottom:12px;",
      metric_card("Samples with at least 1 potential outlier", n_samples_with_outlier),
      metric_card("Candidate outlier values", n_candidate_values)
    ),
    tags$span("Top 10 potential outliers are listed below. The full list of potential outliers (candidate_outliers_*today's_date*.xlsx) is availble for download on tab 4. Export Corrected Data, Plots, Stats, and Report."),
    tags$div(
      style = "border:1px solid #ddd; border-radius:6px; background:white;",
      tags$table(
        style = "border-collapse:collapse; width:100%; font-size:90%;",
        tags$style(htmltools::HTML("
          table, th, td { border: 1px solid #ccc; }
          th, td { padding: 6px 8px; vertical-align: top; }
          th { background:#f7f7f7; text-align:left; }
        ")),
        header, tags$tbody(body_rows)
      )
    )
  )
}

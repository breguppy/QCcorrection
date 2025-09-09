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
ui_filter_info <- function(mv_removed, mv_cutoff) {
  if (length(mv_removed) == 0) {
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
          " metabolites removed based on missing value percentage above ",
          mv_cutoff,
          "%"
        )
      ),
      tags$ul(lapply(mv_removed, tags$li))
    )
  }
  
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


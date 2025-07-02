# Helper function for app.R
library(shiny)
library(bslib)
library(dplyr)
library(shinycssloaders)
library(purrr)

#–– UI snippets ––#
# 1.2 Column-warning text
columnWarningUI <- function(data, selected) {
  warnings <- list()
  
  if (any(selected == "")) {
    warnings[[length(warnings) + 1]] <- tags$span(
      style = "color:darkorange; font-weight:bold;",
      icon("exclamation-triangle"),
      " Please select all four columns."
    )
  } else if (length(unique(selected)) < 4) {
    warnings[[length(warnings) + 1]] <- tags$span(
      style = "color:darkred; font-weight:bold;",
      icon("exclamation-triangle"),
      " Each selected column must be unique."
    )
  }
  
  if (!any(selected == "") && length(unique(selected)) == 4) {
    samp_vec  <- data[[selected[1]]]
    order_vec <- data[[selected[4]]]
    
    if (anyDuplicated(samp_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(
        style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
        icon("exclamation-triangle"),
        " Duplicate sample names detected!"
      )
    }
    
    if (anyDuplicated(order_vec) > 0) {
      warnings[[length(warnings) + 1]] <- tags$span(
        style = "color:darkred; font-weight:bold; display:block; margin-top:5px;",
        icon("exclamation-triangle"),
        " Duplicate order values detected!"
      )
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

# 1.3 Metric card
metric_card <- function(label, value) {
  div(
    style = "background:#f8f9fa; padding:10px; border-radius:8px; flex:1;",
    p(style="font-size:1.4em; font-weight:bold; margin:0;", value),
    h5(style="margin:0;", label)
  )
}

# 1.4 Basic‐info panel
basicInfoUI <- function(df, replacement_counts) {
  metab_cols <- setdiff(names(df), c("sample","batch","class","order"))
  n_metab = length(metab_cols)
  n_missv = sum(is.na(df[,metab_cols]))
  n_qcs   = sum(is.na(df$class))
  n_samp  = sum(!is.na(df$class))
  n_bat   = n_distinct(df$batch)
  n_class = n_distinct(df$class, na.rm=TRUE)
  class_list <- sort(unique(df$class[!is.na(df$class)]))
  
  qc_per_batch <- df %>%
    group_by(batch) %>%
    summarise(qc_in_class = sum(is.na(class)), .groups = "drop")
  
  total_replaced <- sum(replacement_counts$non_numeric_replaced +
                      replacement_counts$zero_replaced)
  
  class_badges <- tags$div(
    style = "display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;",
    lapply(class_list, function(cls) {
      tags$span(
        style = "background-color: #e9ecef; padding: 5px 10px; border-radius: 12px;",
        as.character(cls)
      )
    })
  )
  
  na_table_html <- tags$div(
    style = "margin-top: 15px;",
    tags$h5("Number of QC Samples per Batch"),
    tags$table(
      class = "table table-bordered table-sm",
      tags$thead(
        tags$tr(
          tags$th("Batch"),
          tags$th("QCs in Batch")
        )
      ),
      tags$tbody(
        lapply(1:nrow(qc_per_batch), function(i) {
          tags$tr(
            tags$td(as.character(qc_per_batch$batch[i])),
            tags$td(qc_per_batch$qc_in_class[i])
          )
        })
      )
    )
  )
  
  tagList(
    if (total_replaced > 0) {
      tags$span(style = "color: darkred; font-weight: bold;", 
                paste(total_replaced, "non-numeric or zero metabolite values were removed."))
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
          metric_card("Missing Values", n_missv),
          metric_card("QC Samples", n_qcs),
          metric_card("Samples", n_samp),
          metric_card("Batches", n_bat),
          metric_card("Classes", n_class),
        ),
        
        #bottom left
        tags$div(
          style = "flex: 1; min-width: 250px;",
          tags$h5("Unique Classes"),
          class_badges
        )
        
      ),
      
      # right half
      # QC per batch table column
      tags$div(
        style = "flex: 1; min-width: 250px;",
        tags$h5("Number of QC Samples per Batch"),
        tags$table(
          class = "table table-bordered table-sm",
          tags$thead(
            tags$tr(
              tags$th("Batch"),
              tags$th("QCs in Batch")
            )
          ),
          tags$tbody(
            lapply(1:nrow(qc_per_batch), function(i) {
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

# 1.5 Filter‐info panel
filterInfoUI <- function(removed) {
  if (length(removed)==0) {
    tags$span(style="color:darkgreen;font-weight:bold;",
         "No metabolites removed.")
  } else {
    tagList(
      tags$span(style="color:darkorange;font-weight:bold;",
           paste(length(removed), "metabolites removed based on missing value threshold:")),
      tags$ul(lapply(removed, tags$li))
    )
  }
}

qcMissingValueWarning <- function(df) {
  metab_cols <- setdiff(names(df), c("sample","batch","class","order"))
  qc_idx <- which(is.na(df$class))
  n_missv = sum(is.na(df[qc_idx, metab_cols]))
  
  if (n_missv > 0) {
    tags$span(style = "color:darkred; font-weight:bold;",
      icon("exclamation-triangle"),
      paste(" ", n_missv," values missing from QC samples")
    )
  } else {
    NULL
  }
  
}

correctionInfoUI <- function(imputed_result, imputeM, corMethod) {
  if (corMethod == "QCRFSC"){
    cor_str <- "QC Random Forest Signal Correction (3 seeds x 500 trees)"
  } else if (corMethod == "QCRLSC") {
    cor_str <- "LOESS polynomial fit"
  }
  ui <- list(tags$div(
      style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin-top 15px;",
      metric_card(paste("missing values imputed with", imputed_result$impute_str), imputed_result$n_missv),
      metric_card(cor_str, "Correction Method:")
      )
  )
  do.call(tagList, ui)  
}

postCorFilterInfoUI <- function(filtered_corrected_result) {
  n_removed <- length(filtered_corrected_result$removed_metabolites)
  if(n_removed > 0) {
    ui <- list(
      tags$span(style = "color: darkorange; font-weight: bold;",
                paste(n_removed, "metabolite columns were removed based on QC RSD threshold.")),
      tags$br(),
      tags$span(style = "font-weight: bold;","Removed Columns:"),
      tags$ul(
        lapply(filtered_corrected_result$removed_metabolites, function(name) {
          tags$li(name)
        })
      )
    )
  } else if (n_removed == 0){
    ui <- list(
      tags$span(style = "color: darkgreen; font-weight: bold;",
                "No metabolite columns were removed based on QC RSD threshold.")
    )
  }
  do.call(tagList, ui)
}

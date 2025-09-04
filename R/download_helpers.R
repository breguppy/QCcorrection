#' @keywords internal

#-- Downloads corrected data file.
corrected_file_download <- function(p, d) {
  # Create a new workbook
  wb <- createWorkbook()
  # make column names bold
  bold_style <- createStyle(textDecoration = "Bold")
  # Description text style with orange background.
  style <- createStyle(wrapText = TRUE,
                       valign = "top",
                       fgFill = "#f8cbad")
  
  # Add 0. Raw Data tab
  addWorksheet(wb, "0. Raw Data")
  writeData(
    wb,
    sheet = "0. Raw Data",
    x = d$filtered$df,
    startRow = 3,
    headerStyle = bold_style
  )
  tab0_description <- paste(
    "Tab 0. This tab contains raw peak areas for stable isotope-labeled (13C, deuterium) internal standards (pre-fix of 0-ISTD) and detected metabolites for both pooled quality control (QC) and experimental samples (S).",
    "Data within this tab are formatted for instrument drift correction using the multiply analyzed pooled QC sample.",
    "Batch, class, and order are identifiers for the correction methods, with batch and class set to 1 unless different batches or sample classes are being analyzed within the same run.",
    "The correction settings are shown on the next tab."
  )
  writeData(
    wb,
    sheet = "0. Raw Data",
    x = tab0_description,
    startCol = 1,
    startRow = 1
  )
  mergeCells(wb,
             sheet = "0. Raw Data",
             cols = 1:12,
             rows = 1)
  addStyle(
    wb,
    sheet = "0. Raw Data",
    style = style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  setRowHeights(wb, "0. Raw Data", rows = 1, heights = 60)
  
  # Add 1. Correction Settings tab
  correction_settings_df <- data.frame(
    Settings = c(
      "Sample Column Name",
      "Batch Column Name",
      "Class Column Name",
      "Order Column Name",
      "Missing Value Threshold",
      "QC Missing Value Imputation Method",
      "Sample Missing Value Imputation Method",
      "Correction Method",
      "Remove Imputed Values After Correction?",
      "QC RSD% Threshold",
      "Scaling/Normalization Method",
      "Exclude ISTD in Scaling/Normalization",
      "Keep Corrected QCs"
    ),
    Values = c(
      p$sample_col,
      p$batch_col,
      p$class_col,
      p$order_col,
      paste0(d$filtered$Frule, "%"),
      
      d$imputed$qc_str,
      d$imputed$sam_str,
      d$corrected$str,
      p$remove_imputed,
      paste0(d$filtered_corrected$rsd_cutoff, "%"),
      p$transform,
      p$ex_ISTD,
      p$keep_corrected_qcs
    ),
    stringsAsFactors = FALSE
  )
  addWorksheet(wb, "1. Correction Settings")
  tab1_description <- paste(
    "Tab 1. This tab contains the correction settings along a list of metabolites that were eliminated throughout the process.",
    "The post-QC corrected data are shown on the next tab."
  )
  writeData(
    wb,
    sheet = "1. Correction Settings",
    x = tab1_description,
    startCol = 1,
    startRow = 1
  )
  mergeCells(wb,
             sheet = "1. Correction Settings",
             cols = 1:3,
             rows = 1)
  addStyle(
    wb,
    sheet = "1. Correction Settings",
    style = style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  setRowHeights(wb,
                "1. Correction Settings",
                rows = 1,
                heights = 60)
  writeData(
    wb,
    sheet = "1. Correction Settings",
    x = correction_settings_df,
    startRow = 3,
    startCol = 1,
    headerStyle = bold_style
  )
  # Append withheld from correction columns
  current_col <- 4
  if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold)) {
    withheld_df <- data.frame(
      Columns_Withheld_From_Correction =d$cleaned$withheld_cols,
      stringsAsFactors = FALSE
    )
    writeData(
      wb,
      "1. Correction Settings",
      x = withheld_df,
      startRow = 3,
      startCol = current_col,
      headerStyle = bold_style
    )
    width_vec <- apply(withheld_df, 2, function(x)
      max(nchar(as.character(x)) + 2, na.rm = TRUE))
    width_vec_header <- nchar(colnames(withheld_df)) + 2
    setColWidths(
      wb,
      sheet = "1. Correction Settings",
      cols = current_col,
      widths = pmax(width_vec, width_vec_header)
    )
    current_col <- current_col + 2
  }
  # Append Missing Value Filtered Metabolites
  if (length(d$filtered$mv_removed_cols) > 0) {
    mv_df <- data.frame(
      Missing_Value_Filtered_Metabolites = d$filtered$mv_removed_cols,
      stringsAsFactors = FALSE
    )
    writeData(
      wb,
      "1. Correction Settings",
      x = mv_df,
      startRow = 3,
      startCol = current_col,
      headerStyle = bold_style
    )
    width_vec <- apply(mv_df, 2, function(x)
      max(nchar(as.character(x)) + 2, na.rm = TRUE))
    width_vec_header <- nchar(colnames(mv_df)) + 2
    setColWidths(
      wb,
      sheet = "1. Correction Settings",
      cols = current_col,
      widths = pmax(width_vec, width_vec_header)
    )
    current_col <- current_col + 2
  }
  # Append QC RSD Filtered Metabolites
  if (length(d$filtered_corrected$removed_metabolites) > 0) {
    rsd_df <- data.frame(
      QC_RSD_Filtered_Metabolites = d$filtered_corrected$removed_metabolites,
      stringsAsFactors = FALSE
    )
    writeData(
      wb,
      "1. Correction Settings",
      x = rsd_df,
      startRow = 3,
      startCol = current_col,
      headerStyle = bold_style
    )
    width_vec <- apply(rsd_df, 2, function(x)
      max(nchar(as.character(x)) + 2, na.rm = TRUE))
    width_vec_header <- nchar(colnames(rsd_df)) + 2
    max_width_vec <- pmax(width_vec, width_vec_header)
    setColWidths(wb,
                 sheet = "1. Correction Settings",
                 cols = current_col,
                 widths = max_width_vec)
    current_col <- current_col + 2
  }
  # Append withheld from transformation columns
  if (length(d$transformed$withheld_cols) > 0) {
    ex_trn <- data.frame(
      Exculded_From_Normalization = d$transformed$withheld_cols,
      stringsAsFactors = FALSE
    )
    writeData(
      wb,
      "1. Correction Settings",
      x = ex_trn,
      startRow = 3,
      startCol = current_col,
      headerStyle = bold_style
    )
    width_vec <- apply(ex_trn, 2, function(x)
      max(nchar(as.character(x)) + 2, na.rm = TRUE))
    width_vec_header <- nchar(colnames(ex_trn)) + 2
    max_width_vec <- pmax(width_vec, width_vec_header)
    setColWidths(wb,
                 sheet = "1. Correction Settings",
                 cols = current_col,
                 widths = max_width_vec)
  }
  # Adjust with of correction settings columns.
  width_vec <- apply(correction_settings_df, 2, function(x)
    max(nchar(as.character(x)) + 2, na.rm = TRUE))
  width_vec_header <- nchar(colnames(correction_settings_df)) + 2
  max_width_vec <- pmax(width_vec, width_vec_header)
  setColWidths(wb,
               sheet = "1. Correction Settings",
               cols = 1:2,
               widths = max_width_vec)
  
  # Add 2. Drift Normalized tab
  corrected_df <- d$filtered_corrected$df
  tab2_description <- paste(
    "Tab 2. This tab shows instrument drift corrected values for metabolite levels in experimental samples. Data is corrected using",
    d$corrected$str,
    "FOr each metabolite, this method",
    d$corrected$parameters,
    "This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
    "This corrects for normal instrument drift during the run.",
    "It produces relative metabolite level values in arbitrary units.",
    "For a given metabolite across the entire run, these values average close to 1 for most metabolites.",
    "These data are further normalized on the next tab."
  )
  if (!p$keep_corrected_qcs) {
    corrected_df <- corrected_df[corrected_df$class != "QC", ]
  }
  addWorksheet(wb, "2. Drift Normalized")
  writeData(
    wb,
    sheet = "2. Drift Normalized",
    x = tab2_description,
    startCol = 1,
    startRow = 1
  )
  mergeCells(wb,
             sheet = "2. Drift Normalized",
             cols = 1:16,
             rows = 1)
  addStyle(
    wb,
    sheet = "2. Drift Normalized",
    style = style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  setRowHeights(wb,
                "2. Drift Normalized",
                rows = 1,
                heights = 60)
  writeData(
    wb,
    sheet = "2. Drift Normalized",
    x = corrected_df,
    startRow = 3,
    headerStyle = bold_style
  )
  
  # Add Scaled or Normalized
  transformed_df <- d$transformed$df
  keep_cols <- setdiff(names(transformed_df), d$transformed$withheld_cols)
  if (!p$keep_corrected_qcs) {
    transformed_df <- transformed_df[transformed_df$class != "QC", keep_cols]
  } else {
    transformed_df <- transformed_df[, keep_cols]
  }
  addWorksheet(wb, "3. Scaled or Normalized")
  TRN_description <- paste(
    "Tab 3. This tab shows metabolite level values ratiometrically normalized to total metabolite signal on a per sample basis.",
    "This normalization is done by summing all individual post-QC corrected metabolite level values within a sample (total signal) and then dividing each individual metabolite level value within that sample by the total signal.",
    "This normalization quantifies individual metabolite values across samples based on their proportion to total metabolite load, in arbitrary units, within each individual sample.",
    "These values are displayed on this tab after multiplying by the total number of metabolites present in the sample for easier visualization. Data remain in arbitrary units.",
    "Because arbitary units for a given metabolite quantitatively scale across samples, levels of a given metabolite may be quantitatively compared across samples.",
    "Because unit scaling is different for each metabolite, different metabolites within in a sample cannot be quantitatively compared.",
    "However, because differences in arbitrary unit scaling between samples cancel out by divsion, within-sample metabolite ratios can be quantitatively compared across samples."
  )
  Log2_description <- "Tab 3. This tab shows the log 2 transformed metabolite level values."
  none_description <- "Tab 3. No scaling or normalization method has been applied to the data."
  if (p$transform == "TRN") {
    tab3_description <- TRN_description
  } else if (p$transform == "log2") {
    tab3_description <- Log2_description
  } else {
    tab3_description <- none_description
  }
  writeData(
    wb,
    sheet = "3. Scaled or Normalized",
    x = tab3_description,
    startCol = 1,
    startRow = 1
  )
  mergeCells(wb,
             sheet = "3. Scaled or Normalized",
             cols = 1:22,
             rows = 1)
  addStyle(
    wb,
    sheet = "3. Scaled or Normalized",
    style = style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  setRowHeights(wb,
                "3. Scaled or Normalized",
                rows = 1,
                heights = 80)
  
  writeData(
    wb,
    sheet = "3. Scaled or Normalized",
    x = transformed_df,
    startRow = 3,
    headerStyle = bold_style
  )
  
  # Add 4. Grouped Data Organized
  grouped_data <- group_stats(transformed_df)
  addWorksheet(wb, "4. Grouped Data Organized")
  tab4_description <- paste(
    "Tab 4. This tab shows post-scaled/normalized metabolite level values sorted by group, with group means, standard erorrs (SE), and coefficients of variation (CV) shown.",
    "Because the Metabolomics Core does not perform formal statistical analysis, these statistical analyses are shown for your convenience and quick appraisal of the data.",
    "For publication, data should be analyzed according to the standards of your field, including with the help of a statistician when needed."
  )
  writeData(
    wb,
    sheet = "4. Grouped Data Organized",
    x = tab4_description,
    startCol = 1,
    startRow = 1
  )
  mergeCells(wb,
             sheet = "4. Grouped Data Organized",
             cols = 1:12,
             rows = 1)
  addStyle(
    wb,
    sheet = "4. Grouped Data Organized",
    style = style,
    rows = 1,
    cols = 1,
    gridExpand = TRUE
  )
  setRowHeights(wb,
                "4. Grouped Data Organized",
                rows = 1,
                heights = 60)
  current_row <- 3
  for (group_name in names(grouped_data$group_dfs)) {
    group <- grouped_data$group_dfs[[group_name]]
    group_size <- nrow(group)
    
    writeData(
      wb,
      sheet = "4. Grouped Data Organized",
      x = group,
      startRow = current_row,
      headerStyle = bold_style
    )
    current_row <- current_row + group_size + 1
    group_stats <- grouped_data$group_stats_dfs[[group_name]]
    writeData(
      wb,
      sheet = "4. Grouped Data Organized",
      x = group_stats,
      startRow = current_row,
      startCol = 2,
      headerStyle = bold_style
    )
    current_row <- current_row + 6
  }
  
  # Add. 5. Grouped Data Fold Change (If control group exists)
  if (p$no_control == "FALSE" && p$control_class != "") {
    control_stats <- grouped_data$group_stats_dfs[[p$control_class]]
    fold_change <- fold_changes(transformed_df, control_stats[1, ])
    group_fc_data <- group_stats(fold_change)
    addWorksheet(wb, "5. Grouped Data Fold Change")
    tab5_description <- paste(
      "Tab 5. This tab shows post-ratiometrically normalized metabolite level values expressed in terms of fold change relative to the",
      p$control_class,
      "group mean.",
      "These values have been sorted by group, with group means, standard errors (SE), and coefficients of variation (CV) shown.",
      "Because the Metabolomics Core does not perform formal statistical analysis, these statistical analyses are shown for your convenience and quick appraisal of the data.",
      "For publication, data should be analyzed according to the standards of your field, including with the help of a statistician when needed."
    )
    writeData(
      wb,
      sheet = "5. Grouped Data Fold Change",
      x = tab5_description,
      startCol = 1,
      startRow = 1
    )
    mergeCells(wb,
               sheet = "5. Grouped Data Fold Change",
               cols = 1:12,
               rows = 1)
    addStyle(
      wb,
      sheet = "5. Grouped Data Fold Change",
      style = style,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    setRowHeights(wb,
                  "5. Grouped Data Fold Change",
                  rows = 1,
                  heights = 60)
    current_row <- 3
    for (group_name in names(group_fc_data$group_dfs)) {
      group <- group_fc_data$group_dfs[[group_name]]
      group_size <- nrow(group)
      
      writeData(
        wb,
        sheet = "5. Grouped Data Fold Change",
        x = group,
        startRow = current_row,
        headerStyle = bold_style
      )
      current_row <- current_row + group_size + 1
      group_stats <- group_fc_data$group_stats_dfs[[group_name]]
      writeData(
        wb,
        sheet = "5. Grouped Data Fold Change",
        x = group_stats,
        startRow = current_row,
        startCol = 2,
        headerStyle = bold_style
      )
      current_row <- current_row + 6
    }
    
    # Add. Appendix1. Metaboanalyst Ready as fold change data (if control group exists)
    names(fold_change)[names(fold_change) == "sample"] <- "Sample Name"
    names(fold_change)[names(fold_change) == "class"] <- "Group"
    fold_change$batch <- NULL
    fold_change$order <- NULL
    addWorksheet(wb, "Appendix1. Metaboanalyst Ready")
    writeData(wb, sheet = "Appendix1. Metaboanalyst Ready", x = fold_change)
    setColWidths(wb,
                 sheet = "Appendix1. Metaboanalyst Ready",
                 cols = 1:2,
                 widths = "auto")
    
  } else {
    # Add. Appendix1. Metaboanalyst Ready as normalized data (if there is no control group)
    names(transformed_df)[names(transformed_df) == "sample"] <- "Sample Name"
    names(transformed_df)[names(transformed_df) == "class"] <- "Group"
    transformed_df$batch <- NULL
    transformed_df$order <- NULL
    addWorksheet(wb, "Appendix1. Metaboanalyst Ready")
    writeData(wb, sheet = "Appendix1. Metaboanalyst Ready", x = transformed_df)
    setColWidths(wb,
                 sheet = "Appendix1. Metaboanalyst Ready",
                 cols = 1:2,
                 widths = "auto")
  }
  
  return(wb)
}

#-- Download figure zip
figure_folder_download <- function(p, d) {
  # create temp folder for figures
  tmp_dir <- tempdir()
  fig_dir <- file.path(tmp_dir, "figures")
  if (dir.exists(fig_dir))
    unlink(fig_dir, recursive = TRUE)
  dir.create(fig_dir)
  
  # create RSD figure folder
  rsd_fig_dir <- file.path(fig_dir, "RSD figures")
  if (dir.exists(rsd_fig_dir))
    unlink(rsd_fig_dir, recursive = TRUE)
  dir.create(rsd_fig_dir)
  
  # create PCA figure folder
  pca_fig_dir <- file.path(fig_dir, "PCA plots")
  if (dir.exists((pca_fig_dir)))
    unlink(pca_fig_dir, recursive = TRUE)
  dir.create(pca_fig_dir)
  
  # create metabolite figure folder
  met_fig_dir <- file.path(fig_dir, "metabolite figures")
  if (dir.exists(met_fig_dir))
    unlink(met_fig_dir, recursive = TRUE)
  dir.create(met_fig_dir)
  
  
  rsd_fig <- make_rsd_plot(p, d)
  rsd_path <- file.path(rsd_fig_dir,
                        paste0("rsd_comparison_", p$rsd_cal, ".", p$fig_format))
  if (p$fig_format == "png") {
    ggsave(
      rsd_path,
      plot = rsd_fig,
      width = 7.5,
      height = 4.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  } else if (p$fig_format == "pdf") {
    ggsave(rsd_path,
           plot = rsd_fig,
           width = 7.5,
           height = 4.5,
           units = "in",
           device = grDevices::cairo_pdf)
  }
  
  pca_fig <- make_pca_plot(p, d)
  pca_path <- file.path(pca_fig_dir,
                        paste0("pca_comparison_", p$color_col, ".", p$fig_format))
  if (p$fig_format == "png") {
    ggsave(
      pca_path,
      plot = pca_fig,
      width = 8.333,
      height = 4.417,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  } else if (p$fig_format == "pdf") {
    ggsave(pca_path,
           plot = pca_fig,
           width = 8.333,
           height = 4.417,
           units = "in",
           device = grDevices::cairo_pdf)
  }
  
  # Create metabolite scatter plots
  raw_cols <- setdiff(names(d$filtered$df), c("sample", "batch", "class", "order"))
  cor_cols <- setdiff(names(d$filtered_corrected$df),
                      c("sample", "batch", "class", "order"))
  cols <- intersect(raw_cols, cor_cols)
  n <- length(cols)
  withProgress(message = "Creating figures...", value = 0, {
    for (i in seq_along(cols)) {
      metab <- cols[i]
      fig <- make_met_scatter(d, metab)
      metab <- sanitize_figname(metab)
      path <- file.path(met_fig_dir, paste0(metab, ".", p$fig_format))
      if (p$fig_format == "png") {
        ggsave(
          path,
          plot = fig,
          width = 5,
          height = 5,
          units = "in",
          dpi = 300
        )
      } else if (p$fig_format == "pdf") {
        ggsave(path,
               plot = fig,
               width = 5,
               height = 5,
               units = "in",
               device = grDevices::cairo_pdf)
      }
      incProgress(1 / n, detail = paste("Saved:", metab))
    }
  })
  
  return(list(tmp_dir = tmp_dir, fig_dir = fig_dir))
}

#-- generate correction report
generate_cor_report <- function(p, d, out_dir, template = system.file("app", "report_templates", "report.Rmd", package = "QCcorrection")) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_pdf <- file.path(out_dir, "correction_report.pdf")
  env <- new.env(parent = globalenv())
  
  # create plots
  # get 2 columns for met scatter plots
  if (p$rsd_cal == "met") {
    top2 <- metabolite_rsd(d$filtered$df)  %>%
      select(Metabolite, RSD_NonQC_before = RSD_NonQC) %>%
      inner_join(
        metabolite_rsd(d$filtered_corrected$df) %>%
          select(Metabolite, RSD_NonQC_after = RSD_NonQC),
        by = "Metabolite"
      ) %>%
      mutate(decrease = RSD_NonQC_before - RSD_NonQC_after) %>%
      filter(is.finite(decrease)) %>%
      arrange(desc(decrease)) %>%
      slice_head(n = 2) %>%
      pull(Metabolite)
    
    increased_qc <- metabolite_rsd(d$filtered$df) %>%
      select(Metabolite, RSD_QC_before = RSD_QC) %>%
      inner_join(
        metabolite_rsd(d$filtered_corrected$df) %>%
          select(Metabolite, RSD_QC_after = RSD_QC),
        by = "Metabolite"
      ) %>%
      filter(RSD_QC_after > RSD_QC_before) %>%
      arrange(desc(RSD_QC_after - RSD_QC_before)) %>%
      pull(Metabolite)
  } else {
    top2 <- class_metabolite_rsd(d$filtered$df) %>%
      filter(class != "QC") %>%
      select(Metabolite, RSD_before = RSD) %>%
      inner_join(
        class_metabolite_rsd(d$filtered_corrected$df) %>%
          filter(class != "QC") %>%
          select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) %>%
      mutate(decrease = RSD_before - RSD_after) %>%
      filter(is.finite(decrease)) %>%
      arrange(desc(decrease)) %>%
      distinct(Metabolite, .keep_all = TRUE) %>% 
      slice_head(n = 2) %>%
      pull(Metabolite)
    
    increased_qc <- class_metabolite_rsd(d$filtered$df) %>%
      filter(class == "QC") %>%
      select(Metabolite, RSD_before = RSD) %>%
      inner_join(
        class_metabolite_rsd(d$filtered_corrected$df) %>%
          filter(class == "QC") %>%
          select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) %>%
      filter(RSD_after > RSD_before) %>%
      arrange(desc(RSD_after - RSD_before)) %>%
      pull(Metabolite)
  }
  met1 <- top2[1]
  met2 <- top2[2]
  
  met1_plot <- make_met_scatter(d, met1)
  met2_plot <- make_met_scatter(d, met2)
  rsd_plot <- make_rsd_plot(p, d)
  pca_plot <- make_pca_plot(p, d)
  
  final_data <- list(
    raw_df             = d$cleaned$df,
    replacement_counts = d$cleaned$replacement_counts,
    filtered           = d$filtered,
    filtered_corrected = d$filtered_corrected,
    Frule              = p$Frule,
    post_cor_filter    = p$post_cor_filter,
    rsd_cutoff         = p$rsd_cutoff
  )
  
  # Get descriptions for plots
  descriptions <- list(
    "Withheld Columns" = sprintf(
      "%s%s",
      tagList(
        tags$span(
          style = "font-weight:bold;", "The following columns are non-metabolite columns providing meta-information about the data:"
        ),
        tags$ul(
          lapply(c("sample = Identifies sample name", "batch = Identifies batch (large sample sets are separated into batches)", "class = Identifies sample type", "order = The order in which the sample was injected into the instrument."), function(name) {
            tags$li(name)
          })
        )
      ),
      if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold)) {
        tagList(
          tags$span(
            style = "font-weight:bold;", "The following columns were withheld from correction:"
          ),
          tags$ul(
            lapply(d$cleaned$withheld_cols, function(name) {
              tags$li(name)
            })
          )
        )
      } else {""}
    ),
    "Imputation Description" = sprintf(
      "%s%s%s",
      if (d$imputed$qc_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.<br/>", d$imputed$qc_str)
      } else { "No missing QC values.<br/>"},
      if (d$imputed$sam_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.", d$imputed$sam_str)
      } else { "No missing sample values."},
      if (d$imputed$qc_str == "nothing to impute" && d$imputed$sam_str == "nothing to impute") {
        ""
      } else if (p$remove_imputed == TRUE) {
        "<br/>Imputed values are removed after correction."
      }
      else {""}
    ),
    "Correction Description" = sprintf(
      "Data was corrected using %s. For each metabolite, this method %s This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
      d$corrected$str, d$corrected$parameters
    ),
    "Transformation Description" = sprintf(
      "%s <br/> %s", 
      d$transformed$str,
      if (length(d$transformed$withheld_cols) > 0) {
        tagList(
          tags$span(
            style = "font-weight:bold;", "The following columns are withheld from the transformation:"
          ),
          tags$ul(
            lapply(d$transformed$withheld_cols, function(name) {
              tags$li(name)
            })
          )
        )
      } else {""}
    ),
    "Metabolite Scatter Plots" = sprintf(
      "These plots show a metabolites before and after signal drift correction before any transformation is applied. The two metabolites shown above have the largest decrease in sample variation. The change in variation was determined by calculating relative standard deviation (RSD) for each metabolite %s %s%s A full explanation of RSD is in the next section.",
      if (p$rsd_cal == "class_met") "grouping by sample class." else "",
      if (isTRUE(!p$post_cor_filter)) "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above " else "",
      if (isTRUE(!p$post_cor_filter)) sprintf("%s%%.", p$rsd_cutoff) else ""
    ),
    "RSD Comparison" = sprintf(
      "In these plots, the green indicates RSD decreased after %s, red indicates RSD increased after %s, and gray indicates no change in RSD. For these figures RSD is calculated for each metabolite %s %s%s <br/> %s ",
      if (p$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation",
      if (p$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation",
      if (p$rsd_cal == "class_met") "grouping by sample class." else "",
      if (isTRUE(!p$post_cor_filter)) "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above " else "",
      if (isTRUE(!p$post_cor_filter)) sprintf("%s%%.", p$rsd_cutoff) else "",
      if (length(increased_qc) > 0) {
        sprintf(
          "<br/>The following metabolites increased QC RSD after correction: <br/> %s <br/> More investiagtion is needed to determine if these metabolites should be excluded from the data.",
          paste(increased_qc, collapse = ", ")
        )
      } else {
        ""
      }
    ),
    "PCA Comparison" = sprintf(
      "This PCA plot shows both the raw data and %s data colored by %s.",
      if (p$pca_compare == "filtered_cor_data") "corrected" else "corrected and transformed",
      p$color_col
    )
  )
  
  params <- list(
    title   = "QC Correction Report",
    notes   = p$notes %||% "",
    plots = list(
      "Metabolite Scatter 1" = met1_plot,
      "Metabolite Scatter 2" = met2_plot,
      "RSD Comparison"       = rsd_plot,
      "PCA Comparison"       = pca_plot
    ),
    include = NULL,
    choices = final_data,
    descriptions = descriptions
  )
  
  # Always render HTML first
  html_out <- rmarkdown::render(
    input         = template,
    output_format = "html_document",
    output_file   = file.path(out_dir, "correction_report.html"),
    params        = params,
    envir         = env
  )
  
  # Then print HTML to PDF with Chrome
  if (is.null(pagedown::find_chrome())) {
    stop("Chrome/Chromium not found. Please install Google Chrome or Microsoft Edge.")
  }
  pagedown::chrome_print(input = html_out, output = out_pdf)
  out_pdf
}
# Downloading helpers
library(openxlsx)


#-- Downloads corrected data file.
corrected_file_download <- function(input, rv) {
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
    x = rv$filtered$df,
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
      input$sample_col,
      input$batch_col,
      input$class_col,
      input$order_col,
      paste0(rv$filtered$Frule, "%"),
      
      rv$imputed$qc_str,
      rv$imputed$sam_str,
      rv$corrected$str,
      input$remove_imputed,
      paste0(rv$filtered_corrected$rsd_cutoff, "%"),
      rv$transformed$str,
      input$ex_ISTD,
      input$keep_corrected_qcs
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
  if (isTRUE(input$withhold_cols) && !is.null(input$n_withhold)) {
    withheld_df <- data.frame(
      Columns_Withheld_From_Correction =rv$cleaned$withheld_cols,
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
  if (length(rv$filtered$mv_removed_cols) > 0) {
    mv_df <- data.frame(
      Missing_Value_Filtered_Metabolites = rv$filtered$mv_removed_cols,
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
  if (length(rv$filtered_corrected$removed_metabolites) > 0) {
    rsd_df <- data.frame(
      QC_RSD_Filtered_Metabolites = rv$filtered_corrected$removed_metabolites,
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
  if (length(rv$transformed$withheld_cols) > 0) {
    ex_trn <- data.frame(
      Exculded_From_Normalization = rv$transformed$withheld_cols,
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
  corrected_df <- rv$filtered_corrected$df
  tab2_description <- paste(
    "Tab 2. This tab shows instrument drift corrected values for metabolite levels in experimental samples. Data is corrected using",
    rv$corrected$str,
    "FOr each metabolite, this method",
    rv$corrected$parameters,
    "This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
    "This corrects for normal instrument drift during the run.",
    "It produces relative metabolite level values in arbitrary units.",
    "For a given metabolite across the entire run, these values average close to 1 for most metabolites.",
    "These data are further normalized on the next tab."
  )
  if (!input$keep_corrected_qcs) {
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
  transformed_df <- rv$transformed$df
  keep_cols <- setdiff(names(transformed_df), rv$transformed$withheld_cols)
  if (!input$keep_corrected_qcs) {
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
  if (input$transform == "TRN") {
    tab3_description <- TRN_description
  } else if (input$transform == "log2") {
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
  if (input$no_control == "FALSE" && input$control_class != "") {
    control_stats <- grouped_data$group_stats_dfs[[input$control_class]]
    fold_change <- fold_changes(transformed_df, control_stats[1, ])
    group_fc_data <- group_stats(fold_change)
    addWorksheet(wb, "5. Grouped Data Fold Change")
    tab5_description <- paste(
      "Tab 5. This tab shows post-ratiometrically normalized metabolite level values expressed in terms of fold change relative to the",
      input$control_class,
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


# TODO: Update this function to use the make_*_plots functions.
#-- Download figure zip
figure_folder_download <- function(input, rv) {
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
  
  
  rsd_fig <- make_rsd_plot(input, rv)
  rsd_path <- file.path(rsd_fig_dir,
                        paste0("rsd_comparison_", input$rsd_cal, ".", input$fig_format))
  if (input$fig_format == "png") {
    ggsave(
      rsd_path,
      plot = rsd_fig,
      width = 7.5,
      height = 4.5,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  } else if (input$fig_format == "pdf") {
    ggsave(rsd_path,
           plot = rsd_fig,
           width = 7.5,
           height = 4.5,
           units = "in",
           device = grDevices::cairo_pdf)
  }
  
  pca_fig <- make_pca_plot(input, rv)
  pca_path <- file.path(pca_fig_dir,
                        paste0("pca_comparison_", input$color_col, ".", input$fig_format))
  if (input$fig_format == "png") {
    ggsave(
      pca_path,
      plot = pca_fig,
      width = 8.333,
      height = 4.417,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  } else if (input$fig_format == "pdf") {
    ggsave(pca_path,
           plot = pca_fig,
           width = 8.333,
           height = 4.417,
           units = "in",
           device = grDevices::cairo_pdf)
  }
  
  # Create metabolite scatter plots
  raw_cols <- setdiff(names(rv$filtered$df), c("sample", "batch", "class", "order"))
  cor_cols <- setdiff(names(rv$filtered_corrected$df),
                      c("sample", "batch", "class", "order"))
  cols <- intersect(raw_cols, cor_cols)
  n <- length(cols)
  withProgress(message = "Creating figures...", value = 0, {
    for (i in seq_along(cols)) {
      metab <- cols[i]
      fig <- make_met_scatter(rv, metab)
      metab <- sanitize_figname(metab)
      path <- file.path(met_fig_dir, paste0(metab, ".", input$fig_format))
      if (input$fig_format == "png") {
        ggsave(
          path,
          plot = fig,
          width = 5,
          height = 5,
          units = "in",
          dpi = 300
        )
      } else if (input$fig_format == "pdf") {
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
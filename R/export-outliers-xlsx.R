#' Export outliers excel file
#' @keywords internal
#' @noRd
export_outliers_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  # make column names bold and descriptions with orange backgroud
  bold  <- openxlsx::createStyle(textDecoration = "Bold")
  note  <- openxlsx::createStyle(wrapText = TRUE,
                                 valign = "top",
                                 fgFill = "#f8cbad")
  
  .add_sheet <- function(name) {
    nm <- gsub("[\\[\\]\\*\\?:/\\\\]", "_", name)
    nm <- substr(nm, 1L, 31L)
    openxlsx::addWorksheet(wb, nm)
    nm
  }
  
  if (p$out_data == "filtered_cor_data") {
    df <- d$filtered_corrected$df
  } else {
    df <- d$transformed$df
  }
  res <- detect_qc_aware_outliers(df, group_nonqc_by_class = p$sample_grouping)
  shiny::withProgress(message = "Creating extreme_values_*today's_date*.xlsx...", value = 0, {
    s1 <- .add_sheet("Sample MD")
    txt1 <- paste(
      "Tab 1. This tab shows sample name, group ID, squared Mahalanois distances, Chi-squared quantile cutoff, and True/False extreme value flag.",
      "The squared Mahalanobis distance of the sample is computed in the robust PC score space within each group.",
      "First, each metabolite is standardized with robust centering (median) and scaling (MAD, IQR/1.349, SD, or 1) within each group.",
      "The sandardized metabolites undergoes PCA where we retain PCs to reach at least 80% variance.",
      "Then a robust covariance is computed using the minimum covariance determinant (MCD), Orthogonalized Gnanadesikan\u002DKettenring (OGK), shrinkage, or classical formula depending on sample size and number of PCs retained.",
      "The next tab will show candidate driver metabolites for the flagged extreme value samples."
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:20, rows = 1)
    openxlsx::addStyle(
      wb,
      s1,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s1, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s1,
      x = res$sample_md,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: Sample MD")
    
    s2 <- .add_sheet("Candidates")
    txt2 <- paste(
      "Tab 2. This tab shows sample name, batch, class, order, group ID, metabolite, and robust z-score for all candidate metabolite extreme values.",
      "The robust z-score for each metabolite is computed for each sample using robust centering (median) and scaling (MAD, IQR/1.349, SD, or 1) within each group.",
      "Samples and metabolites are only displayed if the absolute value of the z-score is above 4 (a conservative z-score cutoff to prevent false positives) if QC RSD is stable (<= 20%) or above 5 if QC RSD is borderline (20% < QC RSD <= 30%).",
      "The next tab will show test results confirming potenial extreme values."
    )
    openxlsx::writeData(wb,
                        s2,
                        x = txt2,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s2, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s2,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s2, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s2,
      x = res$candidate_metabolites,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: Candidates")
    
    s3 <- .add_sheet("Confirmations")
    txt3 <- paste(
      "Tab 3. This tab shows sample name, batch, class, order, groupID, metabolite, robust z-score, QC RSD, test method, p value, test strength, and decision on all candidate extreme values from the previous tab.",
      "For each candidate, if the QC RSD is stable (QC RSD <= 20%) the candidate is tested.",
      "If the candidate's QC RSD is borderline (20% < QC RSD <= 30%) and absolute value of the z-score is greater than or equal to 5, the value is tested.",
      "If the QC RSD is unstable (greater than 30%) the candidate is not tested.",
      "Rosner/ESD test is used for sample size n > 25.",
      "Dixon test is used if the candidate is the group's unique min/max and sample size 3 <= n <= 30, otherwise Grubbs test is used to confirm extreme values.",
      "Tied or ineligible cases can still be confirmed when the sample's Mahalanobis distance is flagged (md_only).",
      "The confirmed candidates are POSSIBLE extreme values. Futher investigation should be done before removing the metabolite values."
    )
    openxlsx::writeData(wb,
                        s3,
                        x = txt3,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s3, cols = 1:16, rows = 1)
    openxlsx::addStyle(
      wb,
      s3,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s3, rows = 1, heights = 60)
    openxlsx::writeData(
      wb,
      s3,
      x = res$confirmations,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / 3, detail = "Saved: Confirmations")
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
  })
  return(wb)
}
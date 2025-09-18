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
  
  # TODO: add notes to each sheet
  s1 <- .add_sheet("QC RSD")
  writeData(wb, s1, res$qc_rsd)
  
  s2 <- .add_sheet("Sample MD")
  writeData(wb, s2, res$sample_md)
  
  s3 <- .add_sheet("Candidates")
  writeData(wb, s3, res$candidate_metabolites)
  
  s4 <- .add_sheet("Confirmations")
  writeData(wb, s4, res$confirmations)
  
  wb
}
#' exports missing value report as excel file
#'
#' @keywords internal
#' @noRd
export_corr_xlsx <- function(df1, df2 = NULL, d_type2 = "Corrected", file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
  
  names(df1)[names(df1) == "col1"] <- "Metabolite 1"
  names(df1)[names(df1) == "col2"] <- "Metabolite 2"
  names(df1)[names(df1) == "cor"]  <- "Pearson's r"
  # make column names bold and descriptions with orange background
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
  
  shiny::withProgress(message = "Creating metabolite correlation summary...", value = 0, {
    if (is.null(df2)) N <- 1 else N <- 2
    
    # sheet 1 Metabolite
    s1 <- .add_sheet("Raw Data")
    txt1 <- paste(
      "Tab Raw Data Correlations. Pearson's r values for all metabolite pairs.",
      "If r = -1 or near -1, the pair have a strong negative linear correlation.",
      "If r = 0, there is no correlation and r values near 0 have weak correlations.",
      "If r = 1 or is close to 1, the pair have a strong positive linear correlation." 
    )
    openxlsx::writeData(wb,
                        s1,
                        x = txt1,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:6, rows = 1)
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
      x = df1,
      startRow = 3,
      startCol = 1,
      headerStyle = bold
    )
    shiny::incProgress(1 / N, detail = "Saved: Raw Data Metabolite Correlations.")
    
    if (!is.null(df2)) {
      names(df2)[names(df2) == "col1"] <- "Metabolite 1"
      names(df2)[names(df2) == "col2"] <- "Metabolite 2"
      names(df2)[names(df2) == "cor"]  <- "Pearson's r"
      s2 <- .add_sheet(paste(d_type2, "Data"))
      txt2 <- paste(
        "Tab", d_type2, "Data Correlations. Pearson's r values for all metabolite pairs.",
        "If r = -1 or near -1, the pair have a strong negative linear correlation.",
        "If r = 0, there is no correlation and r values near 0 have weak correlations.",
        "If r = 1 or is close to 1, the pair have a strong positive linear correlation." 
      )
      openxlsx::writeData(wb,
                          s2,
                          x = txt2,
                          startCol = 1,
                          startRow = 1)
      openxlsx::mergeCells(wb, s2, cols = 1:6, rows = 1)
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
        x = df2,
        startRow = 3,
        startCol = 1,
        headerStyle = bold
      )
      shiny::incProgress(1 / N, detail = paste("Saved:", d_type2, "Date Correlations"))
    }
  })
  return(wb)
}
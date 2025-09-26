#' exports corrected data as excel file
#'
#' @keywords internal
#' @noRd
export_stats_xlsx <- function(p, d, file = NULL) {
  .require_pkg("openxlsx", "write Excel workbooks")
  wb <- openxlsx::createWorkbook()
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
  
  df_before <- d$filtered$df
  # Create excel sheet for RSD values
  if (p$rsd_compare == "filtered_cor_data") {
    df_after <- d$filtered_corrected$df
    s2_name <- "Corrected RSD"
  } else {
    df_after <- d$transformed$df
    s2_name <- "Transformed Corrected RSD"
  }
  
  # Compute RSD based on RSD calc
  if (p$rsd_cal == "met") {
    rsdBefore <- metabolite_rsd(df_before)
    rsdAfter <- metabolite_rsd(df_after)
    txt <- "RSD values are computed for each metabolite. For each metabolite, RSD = 100% * (standard deviation) / (mean)."
    # make comparison df
    df_compare <- inner_join(
      rsdBefore %>% rename(RSD_QC_before = RSD_QC, RSD_NonQC_before = RSD_NonQC),
      rsdAfter  %>% rename(RSD_QC_after  = RSD_QC, RSD_NonQC_after  = RSD_NonQC),
      by = "Metabolite"
    ) %>%
      transmute(
        Metabolite,
        delta_RSD_QC    = RSD_QC_after - RSD_QC_before,
        delta_RSD_NonQC = RSD_NonQC_after - RSD_NonQC_before
      )
  } else {
    rsdBefore <- class_metabolite_rsd(df_before)
    rsdAfter <- class_metabolite_rsd(df_after)
    txt <- "RSD values are computed for each metabolite and class. For each metabolite, RSD = 100% * (class standard deviation) / (class mean)."
    df_compare <- inner_join(
      rsdBefore %>% rename(RSD_before = RSD),
      rsdAfter %>% rename(RSD_after  = RSD),
      by = c("class", "Metabolite")
    ) %>%
      transmute(Metabolite, class, delta_RSD = RSD_after - RSD_before)
  }
  
  shiny::withProgress(message = "Creating rsd_stats_*today's_date*.xlsx...", value = 0, {
    # Create excel sheet for RSD values
    s1 <- .add_sheet("Raw RSD")
    openxlsx::writeData(wb,
                        s1,
                        x = rsdBefore,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s1,
                        x = txt,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s1, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s1,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s1, rows = 1, heights = 40)
    
    shiny::incProgress(1 / 3, detail = "Saved: Raw RSD")
    
    s2 <- .add_sheet(s2_name)
    openxlsx::writeData(wb,
                        s2,
                        x = rsdAfter,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s2,
                        x = txt,
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
    openxlsx::setRowHeights(wb, s2, rows = 1, heights = 40)
    shiny::incProgress(1 / 3, detail = "Saved: (Transformed) Corrected RSD")
    # Compare before and after RSD
    s3 <- .add_sheet("RSD Comparison")
    openxlsx::writeData(wb,
                        s3,
                        x = df_compare,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(
      wb,
      s3,
      x = paste(
        txt,
        "The change in RSD (delta = after - before) values will be negative for a decrease in RSD, positive for an increase in RSD, and zero for no change in RSD."
      ),
      startCol = 1,
      startRow = 1
    )
    openxlsx::mergeCells(wb, s3, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s3,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s3, rows = 1, heights = 40)
    shiny::incProgress(1 / 3, detail = "Saved: RSD Comparison")
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
    
  })
  return(wb)
}
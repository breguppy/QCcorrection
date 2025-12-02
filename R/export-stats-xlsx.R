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
    s2_name <- "Corrected"
  } else {
    df_after <- d$transformed$df
    s2_name <- "Transformed Corrected"
  }
  
  # Compute RSD based on RSD calc
  rsdBefore_met <- metabolite_rsd(df_before)
  rsdAfter_met <- metabolite_rsd(df_after)
  txt_met <- "RSD values are computed for each metabolite. For each metabolite, RSD = 100% * (standard deviation) / (mean)."
    # make comparison df
  df_compare_met <- inner_join(
    rsdBefore_met %>% rename(RSD_QC_before = RSD_QC, RSD_NonQC_before = RSD_NonQC),
    rsdAfter_met  %>% rename(RSD_QC_after  = RSD_QC, RSD_NonQC_after  = RSD_NonQC),
    by = "Metabolite"
  ) %>%
    transmute(
      Metabolite,
      delta_RSD_QC    = RSD_QC_after - RSD_QC_before,
      delta_RSD_NonQC = RSD_NonQC_after - RSD_NonQC_before
    )
  rsdBefore_class <- class_metabolite_rsd(df_before)
  rsdAfter_class <- class_metabolite_rsd(df_after)
  txt_class <- "RSD values are computed for each metabolite and class. For each metabolite, RSD = 100% * (class standard deviation) / (class mean)."
  df_compare_class <- inner_join(
    rsdBefore_class %>% rename(RSD_before = RSD),
    rsdAfter_class %>% rename(RSD_after  = RSD),
    by = c("class", "Metabolite")
  ) %>%
    transmute(Metabolite, class, delta_RSD = RSD_after - RSD_before)
  
  shiny::withProgress(message = "Creating rsd_stats_*today's_date*.xlsx...", value = 0, {
    # Create excel sheet for RSD values
    s1 <- .add_sheet("Raw Metabolite RSD")
    openxlsx::writeData(wb,
                        s1,
                        x = rsdBefore_met,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s1,
                        x = txt_met,
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
    
    shiny::incProgress(1 / 6, detail = "Saved: Raw Metabolite RSD")
    
    s2 <- .add_sheet(paste(s2_name, "Metabolite RSD"))
    openxlsx::writeData(wb,
                        s2,
                        x = rsdAfter_met,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s2,
                        x = txt_met,
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
    shiny::incProgress(1 / 6, detail = "Saved: Corrected Metabolite RSD")
    # Compare before and after RSD
    s3 <- .add_sheet("Metabolite RSD Comparison")
    openxlsx::writeData(wb,
                        s3,
                        x = df_compare_met,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(
      wb,
      s3,
      x = paste(
        txt_met,
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
    shiny::incProgress(1 / 6, detail = "Saved: Metabolite RSD Comparison")
    
    
    s4 <- .add_sheet("Raw Class RSD")
    openxlsx::writeData(wb,
                        s4,
                        x = rsdBefore_class,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s4,
                        x = txt_class,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s4, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s4,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s4, rows = 1, heights = 40)
    
    shiny::incProgress(1 / 6, detail = "Saved: Raw Class RSD")
    
    s5 <- .add_sheet(paste(s2_name, "Class RSD"))
    openxlsx::writeData(wb,
                        s5,
                        x = rsdAfter_class,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(wb,
                        s5,
                        x = txt_met,
                        startCol = 1,
                        startRow = 1)
    openxlsx::mergeCells(wb, s5, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s5,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s5, rows = 1, heights = 40)
    shiny::incProgress(1 / 6, detail = "Saved: Corrected Class RSD")
    # Compare before and after RSD
    s6 <- .add_sheet("Class RSD Comparison")
    openxlsx::writeData(wb,
                        s6,
                        x = df_compare_class,
                        startRow = 3,
                        headerStyle = bold)
    openxlsx::writeData(
      wb,
      s6,
      x = paste(
        txt_class,
        "The change in RSD (delta = after - before) values will be negative for a decrease in RSD, positive for an increase in RSD, and zero for no change in RSD."
      ),
      startCol = 1,
      startRow = 1
    )
    openxlsx::mergeCells(wb, s6, cols = 1:12, rows = 1)
    openxlsx::addStyle(
      wb,
      s6,
      style = note,
      rows = 1,
      cols = 1,
      gridExpand = TRUE
    )
    openxlsx::setRowHeights(wb, s6, rows = 1, heights = 40)
    shiny::incProgress(1 / 6, detail = "Saved: Class RSD Comparison")
    
    
    if (!is.null(file)) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      return(normalizePath(file, winslash = "/"))
    }
    
  })
  return(wb)
}
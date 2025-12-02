#' Read input file of raw data
#'
#' @keywords internal
#' @noRd
read_raw_data <- function(file_path) {
  # Get file extension
  file_ext <- tools::file_ext(file_path)
  
  # read data based on file extension
  df <- switch(
    tolower(file_ext),
    "csv" = read.csv(file_path, header = TRUE, check.names = FALSE, na = c("NA", "N/A", "NaN", "")),
    "xls" = read_excel(file_path, na = c("NA", "N/A", "NaN", "")),
    "xlsx" = read_excel(file_path, na = c("NA", "N/A", "NaN", "")),
    stop(
      "Unsupported file type. Please upload a .csv, .xls, or .xlsx file."
    )
  )
  
  return(df)
}
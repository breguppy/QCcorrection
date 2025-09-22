#' Orchestrates correction step
#'
#' @keywords internal
#' @noRd
correct_data <- function(df, metab_cols, corMethod) {
  if (corMethod == "RF") {
    correction_str <- "Random Forest"
    parameters <- "builds 3 models with seeds 42, 31416, 272. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (rf_correction(df, metab_cols, ntree = 500, seed = seed))
      
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "LOESS") {
    correction_str <- "LOESS"
    parameters <- "builds local polynomials of degree 2 that span 0.75 of the total QC values."
    df_corrected <- loess_correction(df, metab_cols)
  } else if (corMethod == "BW_RF") {
    correction_str <- "Batchwise Random Forest"
    parameters <- "3 models with seeds 42, 31416, 272 are built for each metabolite in each batch. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (bw_rf_correction(df, metab_cols, ntree = 500, seed = seed))
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- .median_across_models(df_list, metadata_cols)
  } else if (corMethod == "BW_LOESS") {
    correction_str <- "Batchwise LOESS"
    parameters <- "builds local polynomials of degree 2 that span 0.75 of the total QC values in each batch."
    df_corrected <- bw_loess_correction(df, metab_cols)
  }
  
  return(list(
    df = df_corrected,
    str = correction_str,
    parameters =  parameters
  ))
}
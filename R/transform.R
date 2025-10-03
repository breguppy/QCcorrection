#' Transformation methods for corrected data
#'
#' @keywords internal
#' @noRd
.total_ratio_norm <- function(df, metab_cols) {
  metab_data <- df[, metab_cols, drop = FALSE]
  
  # sum metab_cols values in each row (sample).
  row_sums <- rowSums(metab_data, na.rm = TRUE)
  # compute number of non-missing values in each row (sample).
  non_missing_counts <- rowSums(!is.na(metab_data))
  
  # determine row ratio = (number of non-missing) / (row sum)
  ratios <- non_missing_counts / row_sums
  
  # multiply each row by ratio.
  trn_data <- sweep(metab_data, 1, ratios, FUN = "*")
  
  df[, metab_cols] <- trn_data
  
  return(df)
}


transform_data <- function(df, transform, withheld_cols, ex_ISTD = TRUE) {
  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols <- setdiff(names(df), meta_cols)
  if (ex_ISTD) {
    # Get column names containing ISTD and add them to withheld columns
    istd <- grep("ISTD", metab_cols, value = TRUE)
    itsd <- grep("ITSD", metab_cols, value = TRUE)
    withheld_cols <- c(withheld_cols, istd, itsd)
    metab_cols <- setdiff(metab_cols, withheld_cols)
  }
  
  transformed_df <- df[, c(meta_cols, metab_cols)]
  
  if (transform == "none") {
    transform_str <- "After correction, no scaling or transformations have been applied."
  } else if (transform == "log2") {
    transform_str <- paste(
      "After correction, the log 2 transformation is applied",
      "to metabolite level values. For skewed data, the log transformation helps",
      "normalize the distribution of values."
    )
    transformed_df[metab_cols] <- log(transformed_df[metab_cols], base = 2)
  } else if (transform == "TRN") {
    transform_str <- paste(
      "After correction, metabolite level values are ratiometrically normalized",
      "to total metabolite signal on a per sample basis. This normalization is",
      "done by summing all individual post-QC corrected metabolite level values",
      "within a sample (total signal) and then dividing each individual",
      "metabolite level value within that sample by the total signal. Next, the",
      "individual values are multiplied by the total number of metabolites",
      "present in the sample for easier visualization. This normalization",
      "quantifies individual metabolite values across samples based on their",
      "proportion to total metabolite load, in arbitrary units, within each",
      "individual sample. Because arbitrary units for a given metabolite",
      "quantitatively scale across samples, levels of a given metabolite may be",
      "quantitatively compared across samples. Because unit scaling is different",
      "for each metabolite, different metabolites within in a sample cannot be",
      "quantitatively compared. However, because differences in arbitrary unit",
      "scaling between samples cancel out by division, within-sample metabolite",
      "ratios can be quantitatively compared across samples."
    )
    transformed_df <- .total_ratio_norm(transformed_df, metab_cols)
  }
  return(list(
    df = transformed_df,
    str = transform_str,
    withheld_cols = withheld_cols
  ))
}
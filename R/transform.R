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


transform_data <- function(filtered_corrected, transform, withheld_cols, ex_ISTD = TRUE) {
  df_mv <- filtered_corrected$df_mv
  df_no_mv <- filtered_corrected$df_no_mv
  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols_mv <- setdiff(names(df_mv), meta_cols)
  metab_cols_no_mv <- setdiff(names(df_no_mv), meta_cols)
  
  if (ex_ISTD) {
    # Get column names containing ISTD and add them to withheld columns
    istd_mv <- grepl("^(ISTD|ITSD)", metab_cols_mv, ignore.case = TRUE)
    istd_no_mv <- grepl("^(ISTD|ITSD)", metab_cols_no_mv, ignore.case = TRUE)
    istd_names_mv <- metab_cols_mv[istd_mv]
    istd_names_no_mv <- metab_cols_no_mv[istd_no_mv]
    # Get column names containing ISTD and add them to withheld columns
    withheld_cols_mv <- c(withheld_cols, istd_names_mv)
    metab_cols_mv <- setdiff(metab_cols_mv, withheld_cols_mv)
    withheld_cols_no_mv <- c(withheld_cols, istd_names_no_mv)
    metab_cols_no_mv <- setdiff(metab_cols_no_mv, withheld_cols_no_mv)
  }
  
  transformed_df_mv <- df_mv[, c(meta_cols, metab_cols_mv)]
  transformed_df_no_mv <- df_no_mv[, c(meta_cols, metab_cols_no_mv)]
  
  if (transform == "none") {
    transform_str <- "After correction, no scaling or transformations have been applied."
  } else if (transform == "log2") {
    transform_str <- paste(
      "After correction, the log 2 transformation is applied",
      "to metabolite level values. For skewed data, the log transformation helps",
      "normalize the distribution of values."
    )
    transformed_df_mv[metab_cols_mv] <- log(transformed_df_mv[metab_cols_mv], base = 2)
    transformed_df_no_mv[metab_cols_mv] <- log(transformed_df_no_mv[metab_cols_no_mv], base = 2)
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
    transformed_df_mv <- .total_ratio_norm(transformed_df_mv, metab_cols_mv)
    transformed_df_no_mv <- .total_ratio_norm(transformed_df_no_mv, metab_cols_no_mv)
  }
  return(list(
    df_mv = transformed_df_mv,
    df_no_mv = transformed_df_no_mv,
    str = transform_str,
    withheld_cols_mv = withheld_cols_mv,
    withheld_cols_no_mv = withheld_cols_no_mv
  ))
}
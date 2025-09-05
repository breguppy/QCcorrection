#' Metabolite/value filtering functions
#'
#' @keywords internal
#' @noRd
filter_by_missing <- function(df, metab_cols, Frule) {
  # get metadata columns
  meta_cols <- setdiff(names(df), metab_cols)
  
  # Compute percentage of missing values per metabolite column
  missing_pct <- sapply(df[metab_cols], function(col) {
    mean(is.na(col)) * 100
  })
  
  # Keep only columns with missing percentage <= Frule
  mv_keep_cols <- metab_cols[missing_pct <= Frule]
  
  # list columns removed due to missing value %
  mv_removed_cols <- setdiff(metab_cols, mv_keep_cols)
  
  # filter data by metabolite missing value
  df_filtered <- df[, c(meta_cols, mv_keep_cols)]
  
  return(list(
    df = df_filtered,
    Frule = Frule,
    mv_removed_cols = mv_removed_cols
  ))
}

remove_imputed_from_corrected <- function(raw_df, corrected_df) {
  # Ensure both data frames are the same shape
  if (!all(dim(raw_df) == dim(corrected_df))) {
    stop("Both data frames must have the same dimensions.")
  }
  
  # Return a new corrected_df with values removed where raw_df is NA
  corrected_df[is.na(raw_df)] <- NA
  return(corrected_df)
}

filter_by_qc_rsd <- function(df,
                             rsd_cutoff,
                             metadata_cols = c("sample", "batch", "class", "order")) {
  # Compute RSD
  rsd_df <- metabolite_rsd(df, metadata_cols)
  
  # Identify which metabolites to keep and remove
  keep_metabolites <- rsd_df$Metabolite[is.na(rsd_df$RSD_QC) |
                                          rsd_df$RSD_QC <= rsd_cutoff]
  remove_metabolites <- rsd_df$Metabolite[!is.na(rsd_df$RSD_QC) &
                                            rsd_df$RSD_QC > rsd_cutoff]
  
  # Columns to retain in filtered data
  final_cols <- c(metadata_cols, keep_metabolites)
  filtered_df = df[, final_cols, drop = FALSE]
  
  # Return a list with the filtered data and removed metabolites
  return(
    list(
      df = filtered_df,
      rsd_cutoff = rsd_cutoff,
      removed_metabolites = remove_metabolites
    )
  )
}
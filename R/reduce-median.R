#' Internal function for computing median across models
#'
#' @keywords internal
#' @noRd
.median_across_models <- function(df_list,
                                     metadata_cols = c("Sample", "Batch", "Order", "Treatment")) {
  # Get metabolite column names from the first dataframe
  all_cols <- colnames(df_list[[1]])
  metabolite_cols <- setdiff(all_cols, metadata_cols)
  
  # Add a composite key for alignment
  df_list_keyed <- lapply(df_list, function(df) {
    df$key <- do.call(paste, c(df[metadata_cols], sep = "||"))
    df
  })
  
  # Sort by the composite key
  df_list_sorted <- lapply(df_list_keyed, function(df)
    df[order(df$key), ])
  
  # Sanity check: all keys must align
  key_ref <- df_list_sorted[[1]]$key
  if (!all(sapply(df_list_sorted, function(df)
    all(df$key == key_ref)))) {
    stop("Composite keys are not aligned across all data frames.")
  }
  
  # Extract aligned metabolite matrices
  metabolite_array <- array(sapply(df_list_sorted, function(df)
    as.matrix(df[, metabolite_cols])),
    dim = c(
      nrow(df_list_sorted[[1]]),
      length(metabolite_cols),
      length(df_list_sorted)
    ))
  
  # Compute median along third dimension
  median_matrix <- apply(metabolite_array, c(1, 2), median, na.rm = TRUE)
  
  # Reconstruct result
  median_df <- df_list_sorted[[1]][, metadata_cols]
  median_df <- cbind(median_df, as.data.frame(median_matrix))
  colnames(median_df)[-(1:length(metadata_cols))] <- metabolite_cols
  
  median_df <- median_df[order(median_df$order), ]
  
  return(median_df)
}
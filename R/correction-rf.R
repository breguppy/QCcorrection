#' Random forest correction methods
#'
#' @keywords internal
#' @noRd
rf_correction <- function(df,
                          metab_cols,
                          ntree = 500,
                          seed = NULL) {
  # Set seed (if given a seed)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  df_corrected <- df
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # get index of QCs and samples
    qc_idx <- which(df$class == "QC")
    non_qc_idx <- which(df$class != "QC")
    
    # If there aren't enough QCs, skip correction
    if (length(qc_idx) < 5) {
      warning(paste("Skipping", metab, "- too few QC samples."))
      next
    }
    
    # Build random forest using QC samples
    rf_model <- randomForest::randomForest(x = data.frame(order = df$order[qc_idx]),
                                           y = df[[metab]][qc_idx],
                                           ntree = ntree)
    
    # Predict values for all samples using their order
    predicted <- stats::predict(rf_model, newdata = data.frame(order = df$order))
    
    # Apply correction
    df_corrected[[metab]] <- df[[metab]] / predicted
  }
  
  return(df_corrected)
}

bw_rf_correction <- function(df,
                             metab_cols,
                             ntree = 500,
                             seed = NULL) {
  # set seed (if given a seed)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  df_corrected <- df
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # Get QCs for this batch and metabolite
    for (batch_id in unique(df$batch)) {
      batch_df <- df[df$batch == batch_id, ]
      batch_idx <- which(df$batch == batch_id)
      qc_idx <- which(batch_df$class == "QC")
      
      # Skip if not enough QCs in this batch
      if (length(qc_idx) < 5) {
        warning(sprintf(
          "Skipping %s in batch %s - too few QC samples.",
          metab,
          batch_id
        ))
        next
      }
      
      # Fit random forest using QC samples in this batch
      rf_model <- randomForest::randomForest(
        x = data.frame(order = batch_df$order[qc_idx]),
        y = batch_df[[metab]][qc_idx],
        ntree = ntree
      )
      
      # Predict all values in the batch
      predicted <- stats::predict(rf_model, newdata = data.frame(order = batch_df$order))
      
      # Apply correction
      df_corrected[[metab]][batch_idx] <- df[[metab]][batch_idx] / predicted
    }
  }
  
  return(df_corrected)
}
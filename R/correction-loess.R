#' LOESS correction methods
#'
#' @keywords internal
#' @noRd
loess_correction <- function(df,
                             metab_cols,
                             degree = 2,
                             span = 0.75) {
  df_corrected <- df
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # QC samples
    qc_idx <- which(df$class == "QC")
    
    dat <- data.frame(intensity = df[[metab]][qc_idx], order     = df$order[qc_idx])
    
    # Fit loess (intensity ~ order)
    loess_model <- stats::loess(
      intensity ~ order,
      data   = dat,
      span   = span,
      degree = degree
    )
    
    # Predict values for all samples using their order
    predicted <- stats::predict(loess_model, newdata = data.frame(order = df$order))
    
    # Apply correction
    df_corrected[[metab]] <- df[[metab]] / predicted
  }
  
  return(df_corrected)
}

bw_loess_correction <- function(df,
                                metab_cols,
                                degree = 2,
                                span = 0.75) {
  df_corrected <- df
  
  batches <- unique(df$batch)
  
  for (metab  in metab_cols) {
    # Loess correction for each batch
    for (b in batches) {
      batch_df <- df[df$batch == b, ]
      batch_idx <- which(df$batch == b)
      
      # QC samples in this batch
      qc_idx <- which(batch_df$class == "QC")
      
      if (length(qc_idx) < 5) {
        warning(
          sprintf(
            "Skipping batch '%s' for metabolite '%s': not enough QC samples",
            b,
            metab
          )
        )
        progress_counter <- progress_counter + 1
        setTxtProgressBar(pb, progress_counter)
        next
      }
      
      dat <- data.frame(intensity = batch_df[[metab]][qc_idx], order     = batch_df$order[qc_idx])
      
      # Fit loess (intensity ~ order)
      loess_model <- stats::loess(
        intensity ~ order,
        data   = dat,
        span   = span,
        degree = degree
      )
      
      # Predict on all samples in this batch
      predicted <- stats::predict(loess_model, newdata = data.frame(order = batch_df$order))
      
      # Apply correction only to this batch
      df_corrected[[metab]][batch_idx] <- df[[metab]][batch_idx] / predicted
    }
  }
  return(df_corrected)
}
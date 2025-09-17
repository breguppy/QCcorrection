#' LOESS correction methods
#'
#' @keywords internal
#' @noRd
loess_correction <- function(df,
                             metab_cols,
                             degree = 2,
                             span   = 0.75) {
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") &&
        identical(df$class[nrow(df)], "QC"))) {
    stop("First and last samples must be QCs.")
  }
  
  qcid <- which(df$class == "QC")
  y    <- seq_len(nrow(df))
  df_corrected <- df
  
  impute_half_min_by_class <- function(x, cls) {
    cls <- addNA(as.factor(cls))
    for (lev in levels(cls)) {
      idx  <- which(cls == lev)
      gmin <- suppressWarnings(min(x[idx], na.rm = TRUE))
      if (is.finite(gmin)) {
        x[idx][is.na(x[idx])] <- gmin / 2
      } else {
        allmin <- suppressWarnings(min(x, na.rm = TRUE))
        if (is.finite(allmin)) x[idx][is.na(x[idx])] <- allmin / 2
      }
    }
    x
  }
  
  for (metab in metab_cols) {
    loess_model <- stats::loess(df[[metab]][qcid] ~ qcid,
                                span   = span,
                                degree = degree)
    
    predicted <- stats::predict(loess_model, y)
    # correction
    corrected <- as.numeric(df[[metab]] / predicted)
    
    # clamp: negatives replaced with NA
    corrected[corrected <= 0] <- NA
    
    df_corrected[[metab]] <- corrected
  }
  
  metab_matrix <- as.matrix(df_corrected[metab_cols])
  transposed <- t(metab_matrix)
  knn_result <- impute::impute.knn(transposed,
                                   rowmax = 0.99,
                                   colmax = 0.99,
                                   maxp = 15000)
  imputed_matrix <- t(knn_result$data)
  df_corrected[metab_cols] <- as.data.frame(imputed_matrix)
  
  df_corrected
}

bw_loess_correction <- function(df,
                                metab_cols,
                                span = 0.75,
                                degree = 2,
                                min_qc = 5,
                                clamp_eps = .Machine$double.eps) {
  # helper: class-wise half-min imputation
  impute_half_min_by_class <- function(x, cls) {
    cls <- as.factor(cls)
    for (lev in levels(cls)) {
      idx <- which(cls == lev)
      gmin <- suppressWarnings(min(x[idx], na.rm = TRUE))
      if (is.finite(gmin)) {
        x[idx][is.na(x[idx])] <- gmin / 2
      } else {
        # fallback if an entire class is NA
        allmin <- suppressWarnings(min(x, na.rm = TRUE))
        if (is.finite(allmin)) x[idx][is.na(x[idx])] <- allmin / 2
      }
    }
    x
  }
  
  df_corrected <- df
  batches <- unique(df$batch)
  
  for (metab in metab_cols) {
    for (b in batches) {
      b_idx <- which(df$batch == b)
      sub   <- df[b_idx, , drop = FALSE]
      
      if (!(identical(sub$class[1], "QC") && identical(sub$class[nrow(sub)], "QC"))) {
        stop(sprintf("Batch '%s' must start and end with QC.", b))
      }
      
      qcid <- which(sub$class == "QC")
      if (length(qcid) < min_qc) {
        warning(sprintf("Skipping batch '%s' for '%s': only %d QC rows (< %d).",
                        b, metab, length(qcid), min_qc))
        next
      }
      
      # fit and predict on batch indices
      fit  <- stats::loess(sub[[metab]][qcid] ~ qcid, span = span, degree = degree)
      y    <- seq_len(nrow(sub))
      pred <- stats::predict(fit, y)
      #pred <- pmax(pred, clamp_eps)
      
      # correction
      corrected <- as.numeric(sub[[metab]]) / pred
      
      # remove negatives → NA
      corrected[corrected <= 0] <- NA
      
      # impute NA by class (QC vs samples or any class labels)
      corrected <- impute_half_min_by_class(corrected, as.factor(sub$class))
      
      # write back
      df_corrected[[metab]][b_idx] <- corrected
    }
  }
  df_corrected
}
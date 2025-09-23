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
  
  if (length(qcid) < 2L)
    stop("Need at least 2 QC rows for LOESS.")
  
  y    <- seq_len(nrow(df))
  df_corrected <- df
  
  for (metab in metab_cols) {
    qc_y <- df[[metab]][qcid]
    
    # if no finite variability in QC, skip fit and mark for cleanup
    if (sum(is.finite(qc_y)) < 2L || length(unique(qc_y[is.finite(qc_y)])) < 2L) {
      corrected <- rep(NA_real_, nrow(df))
    } else {
      deg <- min(degree, max(1L, length(qcid) - 1L))
      fit <- stats::loess(qc_y ~ qcid, span = span, degree = deg)
      pred <- suppressWarnings(stats::predict(fit, y))
      pred[!is.finite(pred) | pred <= 0] <- NA_real_
      corrected <- as.numeric(df[[metab]] / pred)
      corrected[!is.finite(corrected) | corrected <= 0] <- NA_real_
    }
    
    df_corrected[[metab]] <- corrected
  }
  
  need_knn <- anyNA(df_corrected[metab_cols])
  if (need_knn) {
    # prefill columns with extreme missingness to avoid knn error
    for (metab in metab_cols) {
      x <- df_corrected[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        min_pos <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        if (is.finite(min_pos)) x[is.na(x)] <- min_pos else x[is.na(x)] <- 0
        df_corrected[[metab]] <- x
      }
    }
    m <- as.matrix(df_corrected[metab_cols])
    m <- t(m)
    kn <- impute::impute.knn(m, rowmax = 1, colmax = 1, maxp = 15000)
    df_corrected[metab_cols] <- as.data.frame(t(kn$data))
  }
  
  # final cleanup: smallest positive per metabolite, else 0; never negative
  df_corrected[metab_cols] <- lapply(df_corrected[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (is.finite(mp)) x[is.na(x)] <- mp else x[is.na(x)] <- 0
    x
  })
  
  df_corrected
}

bw_loess_correction <- function(df,
                                metab_cols,
                                span = 0.75,
                                degree = 2,
                                min_qc = 5) {
  
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
      
      qc_y <- sub[[metab]][qcid]
      if (sum(is.finite(qc_y)) < 2L || length(unique(qc_y[is.finite(qc_y)])) < 2L) {
        corrected <- rep(NA_real_, nrow(sub))
      } else {
        deg <- min(degree, max(1L, length(qcid) - 1L))
        fit <- stats::loess(qc_y ~ qcid, span = span, degree = deg)
        y   <- seq_len(nrow(sub))
        pred <- suppressWarnings(stats::predict(fit, y))
        pred[!is.finite(pred) | pred <= 0] <- NA_real_
        corrected <- as.numeric(sub[[metab]] / pred)
        corrected[!is.finite(corrected) | corrected <= 0] <- NA_real_
      }
      df_corrected[[metab]][b_idx] <- corrected
    }
  }
  
  # KNN if helpful, with prefill to avoid 99%-missing error
  if (anyNA(df_corrected[metab_cols])) {
    for (metab in metab_cols) {
      x <- df_corrected[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        if (is.finite(mp)) x[is.na(x)] <- mp else x[is.na(x)] <- 0
        df_corrected[[metab]] <- x
      }
    }
    m <- t(as.matrix(df_corrected[metab_cols]))
    kn <- impute::impute.knn(m, rowmax = 1, colmax = 1, maxp = 15000)
    df_corrected[metab_cols] <- as.data.frame(t(kn$data))
  }
  
  df_corrected[metab_cols] <- lapply(df_corrected[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    if (is.finite(mp)) x[is.na(x)] <- mp else x[is.na(x)] <- 0
    x
  })
  
  df_corrected
}
#' @keywords internal
#' @noRd
.safe_loess_predict <- function(qcid, qc_y, newx, span, degree) {
  ok <- is.finite(qc_y)
  if (sum(ok) < 2L) {
    # no usable QC points → flat 1s
    return(rep(1, length(newx)))
  }
  qx <- qcid[ok]; qy <- qc_y[ok]
  
  # If there are less than 2 unique QC values: linear interp
  if (length(unique(qy)) < 2L || length(qx) < 2L) {
    return(stats::approx(qx, qy, xout = newx, rule = 2)$y)
  }
  
  # If there are at least 2 unique QC values: 
  # LOESS degree 1 if there are only 2 unique QCs 
  # LOESS degree 2 if there are >= 3 unique QC values
  deg <- min(degree, max(1L, length(qx) - 1L))
  
  # minimumn required span is 3 / number of unique QCs
  spn <- max(span, min(1, 3 / length(qx)))
  
  pred <- tryCatch({
    suppressWarnings({
      fit <- stats::loess(qy ~ qx,
                          span   = spn,
                          degree = deg,
                          control = stats::loess.control(surface = "direct"))
      stats::predict(fit, newx)
    })
  }, error = function(e) NA_real_)
  
  
  # If LOESS fails: linear interpolation.
  if (!is.numeric(pred) || all(!is.finite(pred))) {
    stats::approx(qx, qy, xout = newx, rule = 2)$y
  } else {
    pred
  }
}

loess_correction <- function(df, metab_cols, degree = 2, span = 0.75) {
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC")))
    stop("First and last samples must be QCs.")
  
  qcid <- which(df$class == "QC")
  if (length(qcid) < 2L) stop("Need at least 2 QC rows for LOESS.")
  
  y <- seq_len(nrow(df))
  out <- df
  
  for (metab in metab_cols) {
    qc_y <- df[[metab]][qcid]
    pred <- .safe_loess_predict(qcid, qc_y, newx = y, span = span, degree = degree)
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    corr <- as.numeric(df[[metab]]) / pred
    corr[!is.finite(corr) | corr <= 0] <- NA_real_
    out[[metab]] <- corr
  }
  
  if (anyNA(out[metab_cols])) {
    for (metab in metab_cols) {
      x <- out[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        x[is.na(x)] <- if (is.finite(mp)) mp else 0
        out[[metab]] <- x
      }
    }
    kn <- impute::impute.knn(t(as.matrix(out[metab_cols])), rowmax = 1, colmax = 1, maxp = 15000)
    out[metab_cols] <- as.data.frame(t(kn$data))
  }
  
  out[metab_cols] <- lapply(out[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })
  out
}

#' @keywords internal
#' @noRd
bw_loess_correction <- function(df, metab_cols, span = 0.75, degree = 2, min_qc = 5) {
  out <- df
  for (metab in metab_cols) {
    for (b in unique(df$batch)) {
      b_idx <- which(df$batch == b)
      sub   <- df[b_idx, , drop = FALSE]
      
      if (!(identical(sub$class[1], "QC") && identical(sub$class[nrow(sub)], "QC")))
        stop(sprintf("Batch '%s' must start and end with QC.", b))
      
      qcid <- which(sub$class == "QC")
      if (length(qcid) < min_qc) {
        warning(sprintf("Skipping batch '%s' for '%s': only %d QC rows (< %d).",
                        b, metab, length(qcid), min_qc))
        next
      }
      
      y <- seq_len(nrow(sub))
      pred <- .safe_loess_predict(qcid, sub[[metab]][qcid], newx = y, span = span, degree = degree)
      pred[!is.finite(pred) | pred <= 0] <- NA_real_
      corr <- as.numeric(sub[[metab]]) / pred
      corr[!is.finite(corr) | corr <= 0] <- NA_real_
      out[[metab]][b_idx] <- corr
    }
  }
  
  if (anyNA(out[metab_cols])) {
    for (metab in metab_cols) {
      x <- out[[metab]]
      if (mean(is.na(x)) >= 0.99) {
        mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
        x[is.na(x)] <- if (is.finite(mp)) mp else 0
        out[[metab]] <- x
      }
    }
    kn <- impute::impute.knn(t(as.matrix(out[metab_cols])), rowmax = 1, colmax = 1, maxp = 15000)
    out[metab_cols] <- as.data.frame(t(kn$data))
  }
  
  out[metab_cols] <- lapply(out[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })
  out
}

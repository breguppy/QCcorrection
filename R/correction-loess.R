#' @keywords internal
#' @noRd
.safe_loess_predict <- function(qcid, qc_y, newx, span, degree) {
  ok <- is.finite(qc_y); qx <- qcid[ok]; qy <- qc_y[ok]
  n <- length(qx)
  if (n < 2L) return(rep(1, length(newx)))
  
  # thresholds
  if (n < 8L) {
    return(stats::approx(qx, qy, xout = newx, rule = 2)$y)
  }
  
  deg  <- if (n < 12L) 1L else min(2L, degree)
  spn  <- max(span, min(1, 8 / n))  # ensure enough neighbors
  
  # fit on log scale for multiplicative drift, robust to outliers
  pred <- tryCatch({
    fit <- stats::loess(log(qy) ~ qx, span = spn, degree = deg,
                        family = "symmetric",
                        control = stats::loess.control(surface = "direct"))
    exp(stats::predict(fit, newx))
  }, error = function(e) NA_real_)
  
  if (!is.numeric(pred) || all(!is.finite(pred)))
    stats::approx(qx, qy, xout = newx, rule = 2)$y
  else pred
}

loess_correction <- function(df, metab_cols, degree = 2, span = 0.75, min_qc = 5) {
  df <- df[order(df$order), , drop = FALSE]
  if (!(identical(df$class[1], "QC") && identical(df$class[nrow(df)], "QC")))
    stop("First and last samples must be QCs.")
  
  qcid <- which(df$class == "QC")
  if (length(qcid) < min_qc) stop(sprintf("Need at least %d QC rows for LOESS.", min_qc))
  
  y <- seq_len(nrow(df))
  out <- df
  
  for (metab in metab_cols) {
    qc_y <- df[[metab]][qcid]
    if (all(qc_y <= 0, na.rm = TRUE)) {
      out[[metab]] <- 0
      next
    }
    pred <- .safe_loess_predict(qcid, qc_y, newx = y, span = span, degree = degree)
    pred[!is.finite(pred) | pred <= 0] <- NA_real_
    corr <- as.numeric(df[[metab]]) / pred
    
    # re-anchor scale so QC median in this batch equals 1
    sf <- stats::median(corr[qcid], na.rm = TRUE)
    if (is.finite(sf) && sf > 0) corr <- corr / sf
    
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
    needs_knn <- anyNA(out[metab_cols]) && length(metab_cols) >= 2
    if (needs_knn) {
      kn <- impute::impute.knn(t(as.matrix(out[metab_cols])), rowmax = 1, colmax = 1, maxp = 15000)
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
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
      # new guard
      if (all(sub[[metab]][qcid] <= 0, na.rm = TRUE)) {
        out[[metab]][b_idx] <- 0
        next
      }
      
      y <- seq_len(nrow(sub))
      pred <- .safe_loess_predict(qcid, sub[[metab]][qcid], newx = y, span = span, degree = degree)
      pred[!is.finite(pred) | pred <= 0] <- NA_real_
      corr <- as.numeric(sub[[metab]]) / pred
      
      # re-anchor scale so QC median in this batch equals 1
      sf <- stats::median(corr[qcid], na.rm = TRUE)
      if (is.finite(sf) && sf > 0) corr <- corr / sf
      
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
    needs_knn <- anyNA(out[metab_cols]) && length(metab_cols) >= 2
    if (needs_knn) {
      kn <- impute::impute.knn(
        t(as.matrix(out[metab_cols])),
        rowmax = 1, colmax = 1, maxp = 15000
      )
      out[metab_cols] <- as.data.frame(t(kn$data))
    }
    
  }
  
  out[metab_cols] <- lapply(out[metab_cols], function(x) {
    x[!is.finite(x) | x < 0] <- NA_real_
    mp <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
    x[is.na(x)] <- if (is.finite(mp)) mp else 0
    x
  })
  out
}

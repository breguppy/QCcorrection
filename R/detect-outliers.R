#' Hotelling T^2 on pooled non-QC samples with dual z-scores
#'
#' This function:
#' 1) extracts metabolite columns,
#' 2) log-transforms and scales them using pooled non-QC samples,
#' 3) computes Hotelling-like T^2 (via Mahalanobis distance) on pooled non-QC,
#' 4) flags outlier non-QC samples based on a chi-square cutoff,
#' 5) for outlier samples, computes:
#'    - global z-scores (pooled non-QC),
#'    - class-based z-scores (within each non-QC class),
#' 6) flags metabolite values where both |global z| and |class z|
#'    exceed thresholds.
#'
#' @param df data.frame
#'   Input data with columns for metadata and metabolite intensities.
#' @param meta_cols character
#'   Names of metadata columns to exclude from multivariate analysis.
#'   Default: c("sample", "batch", "class", "order").
#' @param class_col character
#'   Name of the column containing class labels. Default: "class".
#' @param qc_label character
#'   Label in `class_col` indicating QC samples (to be excluded from T^2).
#'   Default: "QC".
#' @param alpha numeric
#'   Significance level for the Hotelling ellipse (default 0.05 => 95% ellipse).
#' @param log_transform logical
#'   Whether to apply log2(x + offset) to metabolite columns. Default TRUE.
#' @param log_offset numeric
#'   Constant added before log-transform to avoid log(0). Default 1.
#' @param z_threshold numeric
#'   Absolute threshold for the global z-score. Default 3.
#' @param class_z_threshold numeric
#'   Absolute threshold for the class-based z-score. Default equal to z_threshold.
#' @param min_complete integer
#'   Minimum number of complete non-QC rows required for covariance estimation.
#'   Default 3.
#' @param drop_constant logical
#'   If TRUE, drop metabolite columns with near-zero variance (based on
#'   complete non-QC rows) before covariance estimation. Default TRUE.
#' @param const_tol numeric
#'   Variance threshold below which a metabolite is considered constant.
#'   Default 1e-12.
#' @param ridge_factor numeric
#'   Factor controlling the ridge size added to the covariance matrix:
#'   ridge = ridge_factor * mean(diag(cov_mat)). Default 1e-6.
#'
#' @return list
#'   A list with components:
#'   - data: original df with added columns:
#'       * T2: squared Mahalanobis distance for non-QC samples (NA for QC or
#'             incomplete rows)
#'       * is_outlier_sample: TRUE if non-QC sample outside (1 - alpha) ellipse
#'       * used_in_fit: TRUE if row used to estimate covariance
#'   - extreme_values: long-format data.frame of flagged metabolite values
#'       * includes metadata, metabolite name, raw/log values,
#'         global and class z-scores, abs values, and T2
#'   - params: list of settings used (alpha, cutoff, thresholds, etc.).
#'
#' @examples
#' # res <- detect_hotelling_nonqc_dual_z(df)
#' @keywords internal
#' @noRd
detect_hotelling_nonqc_dual_z <- function(
    df,
    meta_cols         = c("sample", "batch", "class", "order"),
    class_col         = "class",
    qc_label          = "QC",
    alpha             = 0.05,
    log_transform     = TRUE,
    log_offset        = 1,
    z_threshold       = 3,
    class_z_threshold = z_threshold,
    min_complete      = 3L,
    drop_constant     = TRUE,
    const_tol         = 1e-12,
    ridge_factor      = 1e-6
) {
  # 1. Basic checks -----------------------------------------------------------
  missing_meta <- setdiff(meta_cols, names(df))
  if (length(missing_meta) > 0L) {
    stop("Missing metadata columns in df: ",
         paste(missing_meta, collapse = ", "))
  }
  if (!class_col %in% names(df)) {
    stop("class_col '", class_col, "' not found in df.")
  }
  
  # 2. Identify metabolite columns (numeric, non-metadata) -------------------
  candidate_cols <- setdiff(names(df), meta_cols)
  met_cols <- candidate_cols[vapply(df[candidate_cols], is.numeric, logical(1))]
  if (length(met_cols) == 0L) {
    stop("No numeric metabolite columns found.")
  }
  
  X_raw <- as.matrix(df[, met_cols, drop = FALSE])
  
  # 3. Log-transform ----------------------------------------------------------
  if (log_transform) {
    X_log <- log2(X_raw + log_offset)
  } else {
    X_log <- X_raw
  }
  
  # 4. Define non-QC mask and complete non-QC rows for fit -------------------
  class_vec <- df[[class_col]]
  nonqc_mask <- !is.na(class_vec) & class_vec != qc_label
  
  # Use non-QC + complete metabolite data for covariance/scaling
  complete_nonqc <- nonqc_mask & stats::complete.cases(X_log)
  X_fit <- X_log[complete_nonqc, , drop = FALSE]
  
  if (sum(complete_nonqc) < min_complete) {
    stop("Too few complete non-QC rows to estimate covariance.")
  }
  
  # 5. Optionally drop constant (near-zero variance) metabolites -------------
  if (drop_constant) {
    v <- apply(X_fit, 2L, stats::var, na.rm = TRUE)
    const_mask <- v <= const_tol | is.na(v)
    
    if (any(const_mask)) {
      dropped_mets <- met_cols[const_mask]
      message(
        "Dropping ", length(dropped_mets),
        " metabolite(s) with near-zero variance among non-QC: ",
        paste(dropped_mets, collapse = ", ")
      )
      
      keep_mask <- !const_mask
      met_cols  <- met_cols[keep_mask]
      X_raw     <- X_raw[, keep_mask, drop = FALSE]
      X_log     <- X_log[, keep_mask, drop = FALSE]
      X_fit     <- X_fit[, keep_mask, drop = FALSE]
    }
  }
  
  if (length(met_cols) == 0L) {
    stop("All metabolite columns were dropped as constant; cannot proceed.")
  }
  
  # 6. Scaling based on pooled non-QC fit rows --------------------------------
  X_fit_scaled <- scale(X_fit)
  center_scaled <- attr(X_fit_scaled, "scaled:center")
  scale_scaled  <- attr(X_fit_scaled, "scaled:scale")
  
  scale_all_rows <- function(X_raw, center_raw, scale_raw) {
    sweep(sweep(X_raw, 2L, center_raw, FUN = "-"), 2L, scale_raw, FUN = "/")
  }
  
  X_scaled_all <- scale_all_rows(X_log, center_scaled, scale_scaled)
  p <- ncol(X_scaled_all)
  
  # 7. Covariance (classical) + ridge, on pooled non-QC fit rows -------------
  cov_mat <- stats::cov(X_fit_scaled)
  
  # Ridge regularization
  d <- diag(cov_mat)
  mean_diag <- mean(d)
  if (!is.finite(mean_diag) || mean_diag <= 0) {
    mean_diag <- 1
  }
  ridge <- ridge_factor * mean_diag
  cov_reg <- cov_mat + diag(ridge, nrow(cov_mat))
  
  # 8. Compute T^2 for non-QC complete rows ----------------------------------
  n <- nrow(df)
  T2          <- rep(NA_real_, n)
  used_in_fit <- rep(FALSE, n)
  
  # Non-QC rows with complete data in retained metabolite columns
  complete_nonqc_retained <- nonqc_mask & stats::complete.cases(X_log[, met_cols, drop = FALSE])
  
  if (any(complete_nonqc_retained)) {
    idx_nonqc <- which(complete_nonqc_retained)
    X_nonqc_scaled <- X_scaled_all[idx_nonqc, , drop = FALSE]
    
    T2[idx_nonqc] <- stats::mahalanobis(
      x      = X_nonqc_scaled,
      center = colMeans(X_fit_scaled),
      cov    = cov_reg
    )
    
    # Rows used in covariance fit
    used_in_fit[which(complete_nonqc)] <- TRUE
  }
  
  # 9. Chi-square cutoff and outlier flag ------------------------------------
  cutoff <- stats::qchisq(1 - alpha, df = p)
  is_outlier_sample <- !is.na(T2) & (T2 > cutoff)
  
  # 10. Global z-scores (pooled non-QC) --------------------------------------
  # These are just X_scaled_all for retained metabolites
  Z_global <- X_scaled_all  # same dimension as X_log for retained mets
  
  # 11. Class-based z-scores (within each non-QC class) ----------------------
  Z_class <- matrix(NA_real_, nrow = nrow(df), ncol = length(met_cols))
  colnames(Z_class) <- met_cols
  
  nonqc_classes <- sort(unique(class_vec[nonqc_mask]))
  
  for (cls in nonqc_classes) {
    idx_cls <- which(class_vec == cls)
    
    X_cls <- X_log[idx_cls, , drop = FALSE]
    
    mu_cls <- apply(X_cls, 2L, mean, na.rm = TRUE)
    sd_cls <- apply(X_cls, 2L, stats::sd,   na.rm = TRUE)
    
    zero_sd <- !is.finite(sd_cls) | sd_cls == 0
    sd_cls[zero_sd] <- NA_real_
    
    # z = (x - mu) / sd (per class)
    Z_cls <- sweep(X_cls, 2L, mu_cls, FUN = "-")
    Z_cls <- sweep(Z_cls, 2L, sd_cls, FUN = "/")
    
    Z_class[idx_cls, ] <- Z_cls
  }
  
  # 12. Flag metabolite values where both |z_global| and |z_class| exceed ----
  #     thresholds, restricted to outlier samples.
  outlier_idx <- which(is_outlier_sample)
  
  extreme_values <- NULL
  
  if (length(outlier_idx) > 0L) {
    # Build masks for global + class-based z thresholds
    Zg <- Z_global[outlier_idx, , drop = FALSE]
    Zc <- Z_class[outlier_idx,  , drop = FALSE]
    
    mask <- !is.na(Zg) & !is.na(Zc) &
      (abs(Zg) >= z_threshold) &
      (abs(Zc) >= class_z_threshold)
    
    idx <- which(mask, arr.ind = TRUE)
    
    if (nrow(idx) > 0L) {
      row_ids_global <- outlier_idx[idx[, "row"]]
      col_ids        <- idx[, "col"]
      
      metabolite_names <- met_cols[col_ids]
      
      extreme_values <- data.frame(
        df[row_ids_global, meta_cols, drop = FALSE],
        metabolite        = metabolite_names,
        value_raw         = X_raw[cbind(row_ids_global, col_ids)],
        value_log         = X_log[cbind(row_ids_global, col_ids)],
        z_global          = Z_global[cbind(row_ids_global, col_ids)],
        abs_z_global      = abs(Z_global[cbind(row_ids_global, col_ids)]),
        z_class           = Z_class[cbind(row_ids_global, col_ids)],
        abs_z_class       = abs(Z_class[cbind(row_ids_global, col_ids)]),
        T2                = T2[row_ids_global],
        stringsAsFactors  = FALSE
      )
    } else {
      extreme_values <- data.frame(
        df[0, meta_cols, drop = FALSE],
        metabolite    = character(0),
        value_raw     = numeric(0),
        value_log     = numeric(0),
        z_global      = numeric(0),
        abs_z_global  = numeric(0),
        z_class       = numeric(0),
        abs_z_class   = numeric(0),
        T2            = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  } else {
    extreme_values <- data.frame(
      df[0, meta_cols, drop = FALSE],
      metabolite    = character(0),
      value_raw     = numeric(0),
      value_log     = numeric(0),
      z_global      = numeric(0),
      abs_z_global  = numeric(0),
      z_class       = numeric(0),
      abs_z_class   = numeric(0),
      T2            = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # 13. data frame 
  out_df <- df
  out_df$T2               <- T2
  out_df$is_outlier_sample <- is_outlier_sample
  out_df$used_in_fit      <- used_in_fit
  
  list(
    data           = out_df,
    extreme_values = extreme_values,
    params         = list(
      alpha             = alpha,
      cutoff            = cutoff,
      z_threshold       = z_threshold,
      class_z_threshold = class_z_threshold,
      log_transform     = log_transform,
      log_offset        = log_offset,
      p                 = p,
      class_col         = class_col,
      qc_label          = qc_label,
      min_complete      = min_complete,
      drop_constant     = drop_constant,
      const_tol         = const_tol,
      ridge_factor      = ridge_factor,
      retained_metabolites = met_cols
    )
  )
}

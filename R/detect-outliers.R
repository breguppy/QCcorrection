#' Hotelling T^2 outlier detection with optional class-wise grouping
#' and handling of singular covariance matrices.
#'
#' This function:
#' 1) extracts metabolite columns,
#' 2) log-transforms and globally scales them,
#' 3) computes (robust) covariance and Hotelling-like T^2 via Mahalanobis
#'    distance, either globally or within groups defined by `group_col`,
#'    with a small ridge added to avoid singular covariance matrices,
#' 4) flags samples outside the (1 - alpha) ellipse,
#' 5) within flagged samples, identifies metabolite values with large |z|.
#'
#' @param df data.frame
#'   Input data with metadata columns and metabolite columns.
#' @param meta_cols character
#'   Names of metadata columns to exclude from multivariate analysis.
#'   Default assumes "sample", "batch", "class", "order".
#' @param alpha numeric
#'   Significance level for the Hotelling ellipse (default 0.05 => 95% ellipse).
#' @param log_transform logical
#'   Whether to apply log2(x + offset) to metabolite columns.
#' @param log_offset numeric
#'   Constant added before log-transform to avoid log(0).
#' @param robust logical
#'   If TRUE, use rrcov::CovMcd for robust covariance; otherwise classical cov().
#' @param z_threshold numeric
#'   Absolute z-score threshold for per-metabolite flags within outlier samples.
#' @param group_col character or NULL
#'   Optional column name in `df` used to group samples for T^2 computation,
#'   e.g. "class". If NULL, T^2 is computed globally across all samples.
#' @param min_complete_per_group integer
#'   Minimum number of complete rows (no NA in metabolite columns) required
#'   per group to estimate covariance.
#' @param drop_constant logical
#'   If TRUE, drop metabolite columns with (near) zero variance (based on
#'   globally complete rows) before covariance estimation.
#' @param const_tol numeric
#'   Variance threshold below which a metabolite is considered constant.
#' @param ridge_factor numeric
#'   Factor controlling the ridge size. Ridge is computed as
#'   `ridge_factor * mean(diag(cov_mat))` and added to the diagonal to ensure
#'   invertibility.
#'
#' @return list
#'   A list with components:
#'   - data: original df with added columns:
#'       * T2: squared Mahalanobis distance
#'       * is_outlier_sample: TRUE if outside (1 - alpha) ellipse
#'       * used_in_fit: TRUE if row was used to estimate covariance
#'   - extreme_values: long-format data.frame of flagged metabolite values
#'       * includes metadata, metabolite name, raw/log/scaled values, |z|, T2
#'   - params: list of settings used.
#'
#' @examples
#' # Global:
#' # res <- detect_hotelling_with_metabolite_flags_grouped(df)
#' # By class:
#' # res <- detect_hotelling_with_metabolite_flags_grouped(df, group_col = "class")
#' @keywords internal
#' @noRd
detect_hotelling_with_metabolite_flags_grouped <- function(
    df,
    meta_cols              = c("sample", "batch", "class", "order"),
    alpha                  = 0.05,
    log_transform          = TRUE,
    log_offset             = 1,
    robust                 = TRUE,
    z_threshold            = 3,
    group_col              = NULL,
    min_complete_per_group = 3L,
    drop_constant          = TRUE,
    const_tol              = 1e-12,
    ridge_factor           = 1e-6
) {
  ## 1. Check metadata columns
  missing_meta <- setdiff(meta_cols, names(df))
  if (length(missing_meta) > 0L) {
    stop("Missing metadata columns in df: ", paste(missing_meta, collapse = ", "))
  }
  
  ## 2. Identify metabolite columns: non-metadata numeric columns
  candidate_cols <- setdiff(names(df), meta_cols)
  met_cols <- candidate_cols[vapply(df[candidate_cols], is.numeric, logical(1))]
  if (length(met_cols) == 0L) {
    stop("No numeric metabolite columns found.")
  }
  
  X_raw <- as.matrix(df[, met_cols, drop = FALSE])
  
  ## 3. Log-transform if requested
  if (log_transform) {
    X_log <- log2(X_raw + log_offset)
  } else {
    X_log <- X_raw
  }
  
  ## 4. Global complete cases for scaling / variance checks
  complete_global <- stats::complete.cases(X_log)
  X_complete_global <- X_log[complete_global, , drop = FALSE]
  
  if (nrow(X_complete_global) < 3L) {
    stop("Too few complete rows to estimate global scaling.")
  }
  
  ## 5. Optionally drop constant (near-zero variance) metabolite columns
  if (drop_constant) {
    v <- apply(X_complete_global, 2L, stats::var, na.rm = TRUE)
    const_mask <- v <= const_tol | is.na(v)
    
    if (any(const_mask)) {
      dropped_mets <- met_cols[const_mask]
      message(
        "Dropping ", length(dropped_mets), " metabolite(s) with near-zero variance: ",
        paste(dropped_mets, collapse = ", ")
      )
      
      # Keep only non-constant metabolites
      keep_mask <- !const_mask
      met_cols  <- met_cols[keep_mask]
      X_raw     <- X_raw[, keep_mask, drop = FALSE]
      X_log     <- X_log[, keep_mask, drop = FALSE]
      X_complete_global <- X_complete_global[, keep_mask, drop = FALSE]
    }
  }
  
  if (length(met_cols) == 0L) {
    stop("All metabolite columns were dropped as constant; cannot proceed.")
  }
  
  ## 6. Global scaling (center + sd from complete rows)
  X_scaled_complete <- scale(X_complete_global)
  center_scaled     <- attr(X_scaled_complete, "scaled:center")
  scale_scaled      <- attr(X_scaled_complete, "scaled:scale")
  
  scale_all_rows <- function(X_raw, center_raw, scale_raw) {
    sweep(sweep(X_raw, 2L, center_raw, FUN = "-"), 2L, scale_raw, FUN = "/")
  }
  
  X_scaled_all <- scale_all_rows(X_log, center_scaled, scale_scaled)
  p <- ncol(X_scaled_all)
  
  ## 7. Prepare output vectors
  n <- nrow(df)
  T2          <- rep(NA_real_, n)
  used_in_fit <- rep(FALSE, n)
  
  ## 8. Grouping logic
  if (is.null(group_col)) {
    group_factor <- factor(rep("all", n))
  } else {
    if (!group_col %in% names(df)) {
      stop("group_col '", group_col, "' not found in df.")
    }
    group_factor <- as.factor(df[[group_col]])
  }
  
  group_levels <- levels(group_factor)
  
  if (robust && !requireNamespace("rrcov", quietly = TRUE)) {
    stop("Package 'rrcov' is required for robust = TRUE. Please install it.")
  }
  
  ## 9. Function to add ridge to covariance and ensure invertibility
  regularize_cov <- function(cov_mat, ridge_factor) {
    d <- diag(cov_mat)
    mean_diag <- mean(d)
    if (!is.finite(mean_diag) || mean_diag <= 0) {
      # fallback: use 1 as scale if covariance is degenerate
      mean_diag <- 1
    }
    ridge <- ridge_factor * mean_diag
    cov_mat + diag(ridge, nrow(cov_mat))
  }
  
  ## 10. Loop over groups and compute T^2
  for (g in group_levels) {
    group_idx <- which(group_factor == g)
    if (length(group_idx) == 0L) next
    
    # complete cases within this group (no NA in metabolite columns)
    group_complete_mask <- stats::complete.cases(X_log[group_idx, , drop = FALSE])
    group_complete_idx  <- group_idx[group_complete_mask]
    
    if (length(group_complete_idx) < min_complete_per_group) {
      next
    }
    
    X_group_scaled_complete <- X_scaled_all[group_complete_idx, , drop = FALSE]
    
    # Covariance and center (robust or classical) in scaled space
    if (robust) {
      cov_fit <- rrcov::CovMcd(X_group_scaled_complete)
      center_g <- cov_fit@center
      cov_g    <- cov_fit@cov
    } else {
      center_g <- colMeans(X_group_scaled_complete)
      cov_g    <- stats::cov(X_group_scaled_complete)
    }
    
    # Regularize covariance to avoid singular matrix
    cov_g_reg <- regularize_cov(cov_g, ridge_factor = ridge_factor)
    
    # Mahalanobis T^2 for complete rows in this group
    T2[group_complete_idx] <- stats::mahalanobis(
      x      = X_scaled_all[group_complete_idx, , drop = FALSE],
      center = center_g,
      cov    = cov_g_reg
    )
    
    used_in_fit[group_complete_idx] <- TRUE
  }
  
  ## 11. Chi-square cutoff for (1 - alpha) ellipse (df = p after any drops)
  cutoff <- stats::qchisq(1 - alpha, df = p)
  is_outlier_sample <- !is.na(T2) & (T2 > cutoff)
  
  ## 12. Augmented data frame
  out_df <- df
  out_df$T2               <- T2
  out_df$is_outlier_sample <- is_outlier_sample
  out_df$used_in_fit      <- used_in_fit
  
  ## 13. Per-metabolite extreme values based on global z-scores in outliers
  Z <- X_scaled_all  # globally scaled z-scores (for retained metabolites only)
  
  is_outlier_mat <- matrix(is_outlier_sample,
                           nrow = nrow(Z),
                           ncol = ncol(Z),
                           byrow = FALSE)
  
  mask <- is_outlier_mat & !is.na(Z) & (abs(Z) >= z_threshold)
  idx  <- which(mask, arr.ind = TRUE)
  
  if (nrow(idx) > 0L) {
    row_ids <- idx[, "row"]
    col_ids <- idx[, "col"]
    
    metabolite_names <- met_cols[col_ids]
    
    extreme_values <- data.frame(
      df[row_ids, meta_cols, drop = FALSE],
      metabolite   = metabolite_names,
      value_raw    = X_raw[cbind(row_ids, col_ids)],
      value_log    = if (log_transform) X_log[cbind(row_ids, col_ids)] else NA_real_,
      value_scaled = Z[cbind(row_ids, col_ids)],
      abs_z        = abs(Z[cbind(row_ids, col_ids)]),
      T2           = T2[row_ids],
      stringsAsFactors = FALSE
    )
  } else {
    extreme_values <- data.frame(
      df[0, meta_cols, drop = FALSE],
      metabolite   = character(0),
      value_raw    = numeric(0),
      value_log    = numeric(0),
      value_scaled = numeric(0),
      abs_z        = numeric(0),
      T2           = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  ## 14. Return results
  list(
    data           = out_df,
    extreme_values = extreme_values,
    params         = list(
      alpha                  = alpha,
      cutoff                 = cutoff,
      z_threshold            = z_threshold,
      log_transform          = log_transform,
      log_offset             = log_offset,
      robust                 = robust,
      p                      = p,
      group_col              = group_col,
      min_complete_per_group = min_complete_per_group,
      drop_constant          = drop_constant,
      const_tol              = const_tol,
      ridge_factor           = ridge_factor,
      retained_metabolites   = met_cols
    )
  )
}

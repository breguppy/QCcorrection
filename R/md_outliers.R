#' Mahalanobis outliers before/after, with per-metabolite contributions
#' @keywords internal
#' @noRd
md_outliers_by_group <- function(p, before, after,
                                 color_col = p$color_col %||% "class",
                                 alpha = 0.975,
                                 use_pcs = FALSE,   # TRUE = use PC1:2 space
                                 k_pcs  = 2,
                                 top_m  = 10,
                                 robust = FALSE,    # needs robustbase
                                 ridge  = 1e-6) {
  
  stopifnot(is.list(before), is.list(after), "df" %in% names(before), "df" %in% names(after))
  meta_cols  <- c("sample","batch","class","order")
  metab_cols <- intersect(setdiff(names(before$df), meta_cols),
                          setdiff(names(after$df),  meta_cols))
  if (length(metab_cols) < 2) stop("Need >=2 metabolite columns.")
  
  # Impute 'after' like your PCA
  aft_df <- if (any(is.na(after$df[metab_cols])))
    impute_missing(after$df, metab_cols, p$qcImputeM, p$samImputeM)$df else after$df
  
  datasets <- list(
    list(name = "Before", df = before$df),
    list(name = "After",  df = aft_df)
  )
  
  # helper: safe covariance with optional robustness and ridge
  cov_fit <- function(X) {
    if (robust) {
      if (!requireNamespace("robustbase", quietly = TRUE))
        stop("robustbase not installed; set robust=FALSE or install it.")
      fit <- robustbase::covMcd(X)
      S <- fit$cov; mu <- fit$center
    } else {
      mu <- colMeans(X)
      S  <- stats::cov(X)
    }
    # ridge if near-singular
    eig <- tryCatch(eigen(S, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA)
    if (any(!is.finite(eig)) || any(eig <= 0)) S <- S + diag(ridge, ncol(S))
    list(mu = mu, S = S)
  }
  
  # main loop per dataset
  out <- lapply(datasets, function(dset) {
    df <- dset$df
    grp <- df[[color_col]]
    if (is.null(grp)) stop(sprintf("Column '%s' not found.", color_col))
    
    # analysis matrix: either standardized metabolites (default) or PC scores
    if (use_pcs) {
      pr  <- stats::prcomp(df[, metab_cols], center = TRUE, scale. = TRUE)
      X0  <- pr$x[, seq_len(min(k_pcs, ncol(pr$x))), drop = FALSE]
      colnames(X0) <- paste0("PC", seq_len(ncol(X0)))
      varnames <- colnames(X0) # contributions will be per-PC, not metabolites
      base_mat <- X0           # used for MD and contributions
    } else {
      Z   <- scale(df[, metab_cols], center = TRUE, scale = TRUE)
      varnames <- colnames(Z)
      base_mat <- Z
    }
    
    # per-group MD
    split_idx <- split(seq_len(nrow(df)), grp)
    res_list <- lapply(names(split_idx), function(g) {
      idx <- split_idx[[g]]
      Xg  <- base_mat[idx, , drop = FALSE]
      if (nrow(Xg) < ncol(Xg) + 1) {
        # too few samples for stable covariance; skip with NA
        return(tibble::tibble(
          dataset = dset$name, group = g,
          sample  = df$sample[idx],
          MD2     = NA_real_,
          cutoff  = NA_real_,
          is_outlier = NA,
          top_contributors = replicate(length(idx), list(character()), simplify = FALSE)
        ))
      }
      
      fit <- cov_fit(Xg)
      MD2 <- stats::mahalanobis(Xg, center = fit$mu, cov = fit$S)
      
      cutoff <- stats::qchisq(alpha, df = ncol(Xg))
      flags  <- MD2 > cutoff
      
      # per-variable contributions for flagged rows
      invS <- tryCatch(solve(fit$S), error = function(e) MASS::ginv(fit$S))
      # contributions: c = v * (invS %*% v), elementwise; sum(c) = MD^2
      contrib_list <- vector("list", length(idx))
      if (any(flags) && !use_pcs) {
        V <- sweep(Xg, 2, fit$mu, FUN = "-")  # n x p
        for (k in which(flags)) {
          v  <- as.numeric(V[k, , drop = FALSE])
          cj <- v * as.numeric(invS %*% v)
          ord <- order(abs(cj), decreasing = TRUE, na.last = NA)
          top <- head(varnames[ord], top_m)
          contrib_list[[k]] <- top
        }
      } else if (any(flags) && use_pcs) {
        # if using PCs, report which PCs dominate the MD
        V <- sweep(Xg, 2, fit$mu, FUN = "-")
        for (k in which(flags)) {
          v  <- as.numeric(V[k, , drop = FALSE])
          cj <- v * as.numeric(invS %*% v)
          ord <- order(abs(cj), decreasing = TRUE, na.last = NA)
          contrib_list[[k]] <- head(varnames[ord], min(top_m, length(varnames)))
        }
      } else {
        # no flags or cannot compute contributions
        contrib_list[] <- list(character())
      }
      
      tibble::tibble(
        dataset = dset$name,
        group   = g,
        sample  = df$sample[idx],
        MD2     = as.numeric(MD2),
        cutoff  = cutoff,
        is_outlier = as.logical(flags),
        top_contributors = contrib_list
      )
    })
    
    dplyr::bind_rows(res_list)
  })
  
  dplyr::bind_rows(out)
}

#' @keywords internal
#' @noRd
detect_qc_aware_outliers <- function(
    df,
    z_threshold = 3.5,
    qc_rsd_stable = 0.20,      # 20% RSD or lower = stable
    qc_rsd_unstable = 0.30,    # >30% RSD = unstable, skip confirm
    alpha = 0.05,
    confirm_method = c("dixon","grubbs"),
    min_group_n = 4,           # need >=4 for Dixon/Grubbs to be meaningful
    md_cutoff_quantile = 0.975 # chi-square cutoff quantile
) {
  confirm_method <- match.arg(confirm_method)
  stopifnot(all(c("sample","batch","class","order") %in% names(df)))
  
  # packages (optional but recommended)
  has_rb <- requireNamespace("robustbase", quietly = TRUE)
  has_out <- requireNamespace("outliers",   quietly = TRUE)
  
  meta_cols <- c("sample","batch","class","order")
  met_cols  <- setdiff(names(df), meta_cols)
  
  #--- helper: robust center/scale by class for z-scores
  med <- function(x) stats::median(x, na.rm = TRUE)
  madn <- function(x) stats::mad(x, constant = 1.4826, na.rm = TRUE) # normal-consistent
  
  #--- QC RSD per metabolite (over class == "QC")
  qc_idx <- which(df$class == "QC")
  qc_rsd <- setNames(rep(NA_real_, length(met_cols)), met_cols)
  if (length(qc_idx) >= 2) {
    for (m in met_cols) {
      v <- df[[m]][qc_idx]
      mu <- mean(v, na.rm = TRUE)
      sdv <- stats::sd(v, na.rm = TRUE)
      qc_rsd[m] <- if (isTRUE(mu != 0) && is.finite(mu)) 100 * sdv / abs(mu) else NA_real_
    }
  }
  qc_rsd_tbl <- data.frame(metabolite = met_cols, qc_rsd = qc_rsd, row.names = NULL)
  
  #--- robust Mahalanobis per class (dimension-safe)
  chisq_cut <- function(k, q = md_cutoff_quantile) stats::qchisq(q, df = k)
  
  md_results <- lapply(split(df, df$class), function(dg) {
    X <- as.matrix(dg[, met_cols, drop = FALSE])
    # drop all-NA columns in this group
    keep <- which(colSums(!is.na(X)) >= 2)
    if (length(keep) == 0) {
      return(data.frame(sample = dg$sample, class = dg$class, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    }
    X <- X[, keep, drop = FALSE]
    
    # robust standardize by feature using group med/MAD
    m0  <- apply(X, 2, med)
    s0  <- apply(X, 2, madn)
    s0[s0 == 0 | !is.finite(s0)] <- 1
    Z   <- sweep(sweep(X, 2, m0, "-"), 2, s0, "/")
    Z[!is.finite(Z)] <- NA
    
    # handle p > n via PCA to rank r = min(n-1, non-NA rank)
    n <- nrow(Z); p <- ncol(Z)
    # pairwise-complete covariance for PCA
    S <- stats::cov(Z, use = "pairwise.complete.obs")
    eig <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig)) {
      # fallback: classical MD in available dimensions
      k <- min(n - 1, p)
      return(data.frame(sample = dg$sample, class = dg$class, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    }
    vals <- pmax(eig$values, 0)
    vecs <- eig$vectors
    # keep positive-eigen components only
    pos  <- which(vals > 1e-8)
    if (length(pos) == 0) {
      return(data.frame(sample = dg$sample, class = dg$class, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    }
    # variance-capture rule + rank cap
    cumprop <- cumsum(vals[pos]) / sum(vals[pos])
    k_keep  <- max(1, min(length(pos), which(cumprop >= 0.80)[1]))
    k_keep  <- min(k_keep, n - 1)
    if (is.na(k_keep) || k_keep < 1) k_keep <- min(length(pos), n - 1)
    V  <- vecs[, pos[seq_len(k_keep)], drop = FALSE]
    Tscores <- Z %*% V                      # scores matrix (n x k_keep)
    
    # robust covariance of scores
    if (has_rb && n > (k_keep + 1)) {
      mcd <- tryCatch(robustbase::covMcd(Tscores), error = function(e) NULL)
      if (!is.null(mcd) && is.finite(det(mcd$cov))) {
        mu <- mcd$center; Cv <- mcd$cov
      } else {
        mu <- colMeans(Tscores, na.rm = TRUE); Cv <- stats::cov(Tscores, use = "pairwise.complete.obs")
      }
    } else {
      mu <- colMeans(Tscores, na.rm = TRUE); Cv <- stats::cov(Tscores, use = "pairwise.complete.obs")
    }
    # MD
    md <- tryCatch(stats::mahalanobis(Tscores, center = mu, cov = Cv), error = function(e) rep(NA_real_, nrow(Tscores)))
    cut <- chisq_cut(k_keep)
    data.frame(sample = dg$sample, class = dg$class, md = as.numeric(md), cutoff = cut, flagged = is.finite(md) & md > cut)
  })
  md_tbl <- do.call(rbind, md_results)
  
  #--- robust z-scores within class for driver search
  z_tbl <- do.call(rbind, lapply(split(df, df$class), function(dg) {
    Zg <- as.data.frame(dg[meta_cols], stringsAsFactors = FALSE)
    for (m in met_cols) {
      v  <- dg[[m]]
      z  <- (v - med(v)) / madn(v)
      Zg[[m]] <- z
    }
    Zg
  }))
  
  # candidates only for samples flagged by MD
  flagged_samples <- subset(md_tbl, flagged)$sample
  cand_long <- NULL
  if (length(flagged_samples) > 0) {
    z_long <- tidyr::pivot_longer(z_tbl, cols = all_of(met_cols), names_to = "metabolite", values_to = "z")
    z_long <- merge(z_long, df[meta_cols], by = c("sample","batch","class","order"))
    cand_long <- subset(z_long, sample %in% flagged_samples & is.finite(z) & (z >= z_threshold | z <= -z_threshold))
  } else {
    cand_long <- z_long <- data.frame()
  }
  
  #--- QC-aware confirmation with Dixon/Grubbs
  confirm_rows <- list()
  if (nrow(cand_long) > 0) {
    # join QC RSD
    cand_long <- merge(cand_long, qc_rsd_tbl, by = "metabolite", all.x = TRUE)
    
    by_key <- interaction(cand_long$class, cand_long$metabolite, drop = TRUE)
    split_idx <- split(seq_len(nrow(cand_long)), by_key)
    
    for (idx in split_idx) {
      rows <- cand_long[idx, , drop = FALSE]
      g    <- rows$class[1]
      m    <- rows$metabolite[1]
      grp  <- df[df$class == g, c("sample", m), drop = FALSE]
      vals <- grp[[m]]
      names(vals) <- grp$sample
      
      # skip if group too small or all NA
      valid <- which(is.finite(vals))
      if (length(valid) < min_group_n) {
        for (r in seq_len(nrow(rows))) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_,
                                                            decision = "insufficient_n")
        }
        next
      }
      
      # QC gate
      rsd <- rows$qc_rsd[1]
      gate <- ifelse(!is.na(rsd) & rsd <= 100 * qc_rsd_stable, "stable",
                     ifelse(!is.na(rsd) & rsd > 100 * qc_rsd_unstable, "unstable", "borderline"))
      
      for (r in seq_len(nrow(rows))) {
        s   <- rows$sample[r]
        x   <- vals[valid]
        # identify which observed value corresponds to the candidate sample
        if (!(s %in% names(x))) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_,
                                                            decision = "no_value")
          next
        }
        
        if (gate == "unstable") {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_,
                                                            decision = "skip_unstable_qc")
          next
        }
        
        # Strengthen requirement for borderline QC by raising z-threshold at confirm stage
        z_ok <- abs(rows$z[r]) >= ifelse(gate == "borderline", max(z_threshold, 5), z_threshold)
        if (!z_ok) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_,
                                                            decision = "below_borderline_z")
          next
        }
        
        if (!has_out) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_,
                                                            decision = "outliers_pkg_missing")
          next
        }
        
        # run test focused on the extreme in question
        if (confirm_method == "dixon") {
          # Dixon is sensitive to min/max only; test the side that matches candidate
          xv <- sort(x, na.last = NA)
          cand_val <- vals[s]
          # if candidate equals min or max, test; else fall back to Grubbs
          if (identical(cand_val, xv[1]) || identical(cand_val, xv[length(xv)])) {
            tt <- tryCatch(outliers::dixon.test(x, two.sided = TRUE), error = function(e) NULL)
            if (is.null(tt)) {
              confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                                method = "dixon", p_value = NA_real_, decision = "test_error")
            } else {
              decision <- ifelse(tt$p.value < alpha, "confirm", "retain")
              confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                                method = "dixon", p_value = as.numeric(tt$p.value), decision = decision)
            }
          } else {
            # not at extremes -> use Grubbs
            tt <- tryCatch(outliers::grubbs.test(x, two.sided = TRUE), error = function(e) NULL)
            if (is.null(tt)) {
              confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                                method = "grubbs", p_value = NA_real_, decision = "test_error")
            } else {
              decision <- ifelse(tt$p.value < alpha, "confirm", "retain")
              confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                                method = "grubbs", p_value = as.numeric(tt$p.value), decision = decision)
            }
          }
        } else {
          # Grubbs two-sided
          tt <- tryCatch(outliers::grubbs.test(x, two.sided = TRUE), error = function(e) NULL)
          if (is.null(tt)) {
            confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                              method = "grubbs", p_value = NA_real_, decision = "test_error")
          } else {
            decision <- ifelse(stats::pvalue(tt) < alpha, "confirm", "retain")
            confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                              method = "grubbs", p_value = as.numeric(stats::pvalue(tt)), decision = decision)
          }
        }
      }
    }
  }
  confirm_tbl <- if (length(confirm_rows)) do.call(rbind, confirm_rows) else
    data.frame(sample=character(), batch=character(), class=character(), order=numeric(),
               metabolite=character(), z=numeric(), qc_rsd=numeric(),
               method=character(), p_value=numeric(), decision=character(), row.names = NULL)
  
  # outputs
  list(
    qc_rsd = qc_rsd_tbl,
    sample_md = md_tbl,
    candidate_metabolites = cand_long,
    confirmations = confirm_tbl
  )
}

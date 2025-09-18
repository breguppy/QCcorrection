#' @keywords internal
#' @noRd
detect_qc_aware_outliers <- function(
    df,
    group_nonqc_by_class,
    z_threshold = 3.5,
    qc_rsd_stable = 0.20,
    qc_rsd_unstable = 0.30,
    alpha = 0.05,
    confirm_method = c("dixon","grubbs"),
    min_group_n = 3,
    md_cutoff_quantile = 0.975
) {
  confirm_method <- match.arg(confirm_method)
  stopifnot(all(c("sample","batch","class","order") %in% names(df)))
  
  has_rb  <- requireNamespace("robustbase", quietly = TRUE)
  has_out <- requireNamespace("outliers",   quietly = TRUE)
  
  meta_cols <- c("sample","batch","class","order")
  met_cols  <- setdiff(names(df), meta_cols)
  met_cols  <- met_cols[sapply(df[met_cols], is.numeric)]  # ensure numeric only
  
  # ---- NEW: auto-fallback if only one non-QC class level
  nonqc_idx <- which(df$class != "QC")
  n_levels  <- length(unique(df$class[nonqc_idx]))
  effective_grouping <- if (isTRUE(group_nonqc_by_class) && n_levels > 1) "by_class" else "all"
  
  df2 <- df
  df2$group_id <- if (effective_grouping == "by_class") {
    ifelse(df2$class == "QC", "QC", as.character(df2$class))
  } else {
    ifelse(df2$class == "QC", "QC", "all")
  }
  
  # helpers
  med  <- function(x) stats::median(x, na.rm = TRUE)
  # MAD fallback: MAD -> IQR/1.349 -> SD -> 1
  rob_scale <- function(v) {
    d <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(d) || d == 0) d <- stats::IQR(v, na.rm = TRUE) / 1.349
    if (!is.finite(d) || d == 0) d <- stats::sd(v,  na.rm = TRUE)
    if (!is.finite(d) || d == 0) d <- 1
    d
  }
  chisq_cut <- function(k, q = md_cutoff_quantile) stats::qchisq(q, df = k)
  
  # QC RSD per metabolite
  qc_idx <- which(df2$class == "QC")
  qc_rsd <- setNames(rep(NA_real_, length(met_cols)), met_cols)
  if (length(qc_idx) >= 2) {
    for (m in met_cols) {
      v <- df2[[m]][qc_idx]; mu <- mean(v, na.rm = TRUE); sdv <- stats::sd(v, na.rm = TRUE)
      qc_rsd[m] <- if (isTRUE(mu != 0) && is.finite(mu)) 100 * sdv / abs(mu) else NA_real_
    }
  }
  qc_rsd_tbl <- data.frame(metabolite = met_cols, qc_rsd = qc_rsd, row.names = NULL)
  
  # Robust MD per group_id (PCA scores + robust cov)
  md_results <- lapply(split(df2, df2$group_id), function(dg) {
    X <- as.matrix(dg[, met_cols, drop = FALSE])
    keep <- which(colSums(!is.na(X)) >= 2)
    if (length(keep) == 0) {
      return(data.frame(sample = dg$sample, group_id = dg$group_id, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    }
    X <- X[, keep, drop = FALSE]
    m0 <- apply(X, 2, med)
    s0 <- apply(X, 2, rob_scale)
    Z  <- sweep(sweep(X, 2, m0, "-"), 2, s0, "/"); Z[!is.finite(Z)] <- NA
    
    n <- nrow(Z); S <- stats::cov(Z, use = "pairwise.complete.obs")
    eig <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig)) return(data.frame(sample = dg$sample, group_id = dg$group_id, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    vals <- pmax(eig$values, 0); vecs <- eig$vectors; pos <- which(vals > 1e-8)
    if (!length(pos)) return(data.frame(sample = dg$sample, group_id = dg$group_id, md = NA_real_, cutoff = NA_real_, flagged = FALSE))
    cumprop <- cumsum(vals[pos]) / sum(vals[pos])
    k_keep  <- max(1, min(length(pos), which(cumprop >= 0.80)[1])); k_keep <- min(k_keep, n - 1); if (is.na(k_keep) || k_keep < 1) k_keep <- min(length(pos), n - 1)
    V <- vecs[, pos[seq_len(k_keep)], drop = FALSE]; Ts <- Z %*% V
    
    if (has_rb && n > (k_keep + 1)) {
      mcd <- tryCatch(robustbase::covMcd(Ts), error = function(e) NULL)
      if (!is.null(mcd) && is.finite(det(mcd$cov))) { mu <- mcd$center; Cv <- mcd$cov } else { mu <- colMeans(Ts, na.rm=TRUE); Cv <- stats::cov(Ts, use="pairwise.complete.obs") }
    } else {
      mu <- colMeans(Ts, na.rm=TRUE); Cv <- stats::cov(Ts, use="pairwise.complete.obs")
    }
    md <- tryCatch(stats::mahalanobis(Ts, center = mu, cov = Cv), error = function(e) rep(NA_real_, nrow(Ts)))
    cut <- chisq_cut(k_keep)
    data.frame(sample = dg$sample, group_id = dg$group_id, md = as.numeric(md), cutoff = cut, flagged = is.finite(md) & md > cut)
  })
  md_tbl <- do.call(rbind, md_results)
  
  # Robust z by group_id
  z_tbl <- do.call(rbind, lapply(split(df2, df2$group_id), function(dg) {
    Zg <- as.data.frame(dg[, c(meta_cols, "group_id")], stringsAsFactors = FALSE)
    for (m in met_cols) {
      v   <- dg[[m]]
      mu  <- stats::median(v, na.rm = TRUE)
      den <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
      if (!is.finite(den) || den == 0) {
        den <- stats::IQR(v, na.rm = TRUE) / 1.349
      }
      if (!is.finite(den) || den == 0) {
        den <- stats::sd(v, na.rm = TRUE)
      }
      if (!is.finite(den) || den == 0) den <- 1
      Zg[[m]] <- (v - mu) / den
    }
    Zg
  }))
  
  # Candidates: flagged non-QC samples only
  flagged_samples <- subset(md_tbl, flagged & group_id != "QC")$sample
  cand_long <- data.frame()
  
  if (length(flagged_samples) > 0) {
    # z_tbl already has: sample,batch,class,order,group_id + metabolites
    z_long <- tidyr::pivot_longer(
      z_tbl,
      cols = all_of(met_cols),
      names_to = "metabolite",
      values_to = "z"
    )
    cand_long <- subset(
      z_long,
      group_id != "QC" &
        sample %in% flagged_samples &
        is.finite(z) &
        (z >= z_threshold | z <= -z_threshold)
    )
  }
  
  # Confirmation per (group_id, metabolite)
  confirm_rows <- list()
  if (nrow(cand_long) > 0) {
    cand_long <- merge(cand_long, qc_rsd_tbl, by = "metabolite", all.x = TRUE)
    
    by_key <- interaction(cand_long$group_id, cand_long$metabolite, drop = TRUE)
    split_idx <- split(seq_len(nrow(cand_long)), by_key)
    
    for (idx in split_idx) {
      rows <- cand_long[idx, , drop = FALSE]
      gid  <- rows$group_id[1]
      met  <- rows$metabolite[1]
      
      grp  <- df2[df2$group_id == gid & df2$class != "QC", c("sample", met), drop = FALSE]
      vals <- grp[[met]]; names(vals) <- grp$sample
      valid <- which(is.finite(vals))
      if (length(valid) < min_group_n) {
        for (r in seq_len(nrow(rows))) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_, decision = "insufficient_n")
        }
        next
      }
      
      rsd <- rows$qc_rsd[1]
      gate <- ifelse(!is.na(rsd) & rsd <= 100 * qc_rsd_stable, "stable",
                     ifelse(!is.na(rsd) & rsd > 100 * qc_rsd_unstable, "unstable", "borderline"))
      
      for (r in seq_len(nrow(rows))) {
        s <- rows$sample[r]
        x <- vals[valid]
        if (!(s %in% names(x))) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_, decision = "no_value")
          next
        }
        if (gate == "unstable") {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_, decision = "skip_unstable_qc")
          next
        }
        z_ok <- abs(rows$z[r]) >= ifelse(gate == "borderline", max(z_threshold, 5), z_threshold)
        if (!z_ok) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_, decision = "below_borderline_z")
          next
        }
        if (!has_out) {
          confirm_rows[[length(confirm_rows) + 1]] <- cbind(rows[r,],
                                                            method = NA_character_, p_value = NA_real_, decision = "outliers_pkg_missing")
          next
        }
        
        # choose test
        if (confirm_method == "dixon") {
          xv <- sort(x, na.last = NA)
          cand_val <- vals[s]
          at_extreme <- (which.min(x) == which(names(x) == s)) || (which.max(x) == which(names(x) == s))
          if (at_extreme) {
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
      }
    }
  }
  
  confirm_tbl <- if (length(confirm_rows)) do.call(rbind, confirm_rows) else
    data.frame(sample=character(), batch=character(), class=character(), order=numeric(),
               group_id=character(), metabolite=character(), z=numeric(), qc_rsd=numeric(),
               method=character(), p_value=numeric(), decision=character(), row.names = NULL)
  
  list(
    qc_rsd = qc_rsd_tbl,
    sample_md = md_tbl,
    candidate_metabolites = cand_long,
    confirmations = confirm_tbl,
    params = list(group_nonqc_by_class = group_nonqc_by_class,
                  z_threshold = z_threshold,
                  alpha = alpha)
  )
}

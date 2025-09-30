#' @keywords internal
#' @noRd
detect_qc_aware_outliers <- function(df,
                                     group_nonqc_by_class,
                                     z_threshold = 4,
                                     qc_rsd_stable = 0.20,
                                     qc_rsd_unstable = 0.30,
                                     alpha = 0.05,
                                     confirm_method = c("dixon", "grubbs"),
                                     min_group_n = 3,
                                     md_cutoff_quantile = 0.95) {
  confirm_method <- match.arg(confirm_method)
  stopifnot(all(c("sample", "batch", "class", "order") %in% names(df)))
  
  has_rb  <- requireNamespace("robustbase", quietly = TRUE)
  has_out <- requireNamespace("outliers", quietly = TRUE)
  has_env <- requireNamespace("EnvStats", quietly = TRUE)
  
  meta_cols <- c("sample", "batch", "class", "order")
  met_cols  <- setdiff(names(df), meta_cols)
  met_cols  <- met_cols[sapply(df[met_cols], is.numeric)]
  
  # if only one non-QC class exists, group all samples into all group
  nonqc_idx <- which(df$class != "QC")
  n_levels  <- length(unique(df$class[nonqc_idx]))
  effective_grouping <- if (isTRUE(group_nonqc_by_class) &&
                            n_levels > 1)
    "by_class"
  else
    "all"
  
  df2 <- df
  df2$group_id <- if (effective_grouping == "by_class") {
    ifelse(df2$class == "QC", "QC", as.character(df2$class))
  } else {
    ifelse(df2$class == "QC", "QC", "all")
  }
  
  # helpers
  med  <- function(x)
    stats::median(x, na.rm = TRUE)
  # robust scaling: MAD or IQR/1.349 or SD or 1
  rob_scale <- function(v) {
    d <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(d) ||
        d == 0)
      d <- stats::IQR(v, na.rm = TRUE) / 1.349
    if (!is.finite(d) || d == 0)
      d <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(d) || d == 0)
      d <- 1
    d
  }
  
  # QC RSD per metabolite
  qc_idx <- which(df2$class == "QC")
  qc_rsd <- setNames(rep(NA_real_, length(met_cols)), met_cols)
  if (length(qc_idx) >= 2) {
    for (m in met_cols) {
      v <- df2[[m]][qc_idx]
      mu <- mean(v, na.rm = TRUE)
      sdv <- stats::sd(v, na.rm = TRUE)
      qc_rsd[m] <- if (isTRUE(mu != 0) &&
                       is.finite(mu))
        100 * sdv / abs(mu)
      else
        NA_real_
    }
  }
  qc_rsd_tbl <- data.frame(metabolite = met_cols,
                           qc_rsd = qc_rsd,
                           row.names = NULL)
  
  # Robust MD per group (PCA scores + robust cov)
  md_results <- lapply(split(df2, df2$group_id), function(dg) {
    X <- as.matrix(dg[, met_cols, drop = FALSE])
    keep <- which(colSums(!is.na(X)) >= 2)
    if (length(keep) == 0) {
      return(
        data.frame(
          sample = dg$sample,
          group_id = dg$group_id,
          md = NA_real_,
          cutoff = NA_real_,
          flagged = FALSE
        )
      )
    }
    X <- X[, keep, drop = FALSE]
    m0 <- apply(X, 2, med)
    s0 <- apply(X, 2, rob_scale)
    Z  <- sweep(sweep(X, 2, m0, "-"), 2, s0, "/")
    Z[!is.finite(Z)] <- NA
    
    n <- nrow(Z)
    S <- stats::cov(Z, use = "pairwise.complete.obs")
    eig <- tryCatch(
      eigen(S, symmetric = TRUE),
      error = function(e)
        NULL
    )
    if (is.null(eig))
      return(
        data.frame(
          sample = dg$sample,
          group_id = dg$group_id,
          md = NA_real_,
          cutoff = NA_real_,
          flagged = FALSE
        )
      )
    vals <- pmax(eig$values, 0)
    vecs <- eig$vectors
    pos <- which(vals > 1e-8)
    if (!length(pos))
      return(
        data.frame(
          sample = dg$sample,
          group_id = dg$group_id,
          md = NA_real_,
          cutoff = NA_real_,
          flagged = FALSE
        )
      )
    cumprop <- cumsum(vals[pos]) / sum(vals[pos])
    k_keep  <- max(1, min(length(pos), which(cumprop >= 0.80)[1]))
    k_keep <- min(k_keep, n - 1)
    if (is.na(k_keep) || k_keep < 1)
      k_keep <- min(length(pos), n - 1)
    V <- vecs[, pos[seq_len(k_keep)], drop = FALSE]
    Ts <- Z %*% V
    p <- ncol(Ts)
    
    # further cap PCs if sample size is too small for robust MCD
    if (p > 0 && n < 2 * p) {
      keep_p <- max(1, min(p, floor((n - 1) / 2)))
      if (keep_p < p) {
        Ts <- Ts[, seq_len(keep_p), drop = FALSE]
        p  <- keep_p
      }
    }
    
    # choose covariance estimator
    if (p == 1) {
      mu <- colMeans(Ts, na.rm = TRUE)
      Cv <- matrix(stats::var(as.numeric(Ts), na.rm = TRUE), 1, 1)
    } else if (has_rb && n >= 2 * p) {
      # safe to use MCD only when n >= 2p
      mcd <- tryCatch(
        robustbase::covMcd(Ts),
        error = function(e)
          NULL
      )
      if (!is.null(mcd) &&
          is.matrix(mcd$cov) && all(is.finite(mcd$cov))) {
        mu <- mcd$center
        Cv <- mcd$cov
      } else {
        mu <- colMeans(Ts, na.rm = TRUE)
        Cv <- stats::cov(Ts, use = "pairwise.complete.obs")
      }
    } else if (requireNamespace("robustbase", quietly = TRUE)) {
      # OGK is robust and works better when n << p
      ogk <- tryCatch(
        robustbase::covOGK(Ts),
        error = function(e)
          NULL
      )
      if (!is.null(ogk) &&
          is.matrix(ogk$cov) && all(is.finite(ogk$cov))) {
        mu <- ogk$center
        Cv <- ogk$cov
      } else {
        mu <- colMeans(Ts, na.rm = TRUE)
        Cv <- stats::cov(Ts, use = "pairwise.complete.obs")
      }
    } else if (requireNamespace("corpcor", quietly = TRUE)) {
      # shrinkage covariance is well-conditioned for small n
      mu <- colMeans(Ts, na.rm = TRUE)
      Cv <- corpcor::cov.shrink(Ts)
    } else {
      # classical covariance as last resort
      mu <- colMeans(Ts, na.rm = TRUE)
      Cv <- stats::cov(Ts, use = "pairwise.complete.obs")
    }
    
    # ensure positive-definite covariance
    ok <- tryCatch({
      ev <- eigen(Cv, symmetric = TRUE, only.values = TRUE)$values
      all(is.finite(ev)) && min(ev, na.rm = TRUE) > 1e-8
    }, error = function(e)
      FALSE)
    if (!ok) {
      eps <- 1e-6
      Cv  <- Cv + diag(eps, ncol(Ts))
    }
    
    md <- tryCatch(
      stats::mahalanobis(Ts, center = mu, cov = Cv),
      error = function(e)
        rep(NA_real_, nrow(Ts))
    )
    
    cut <- stats::qchisq(md_cutoff_quantile, df = ncol(Ts))
    data.frame(
      sample = dg$sample,
      group_id = dg$group_id,
      md = as.numeric(md),
      cutoff = cut,
      flagged = is.finite(md) & md > cut
    )
  })
  md_tbl <- do.call(rbind, md_results)
  
  # Robust z by group using MAD, IQR / 1.349, SD, or 1 as denominator.
  z_tbl <- do.call(rbind, lapply(split(df2, df2$group_id), function(dg) {
    Zg <- as.data.frame(dg[, c(meta_cols, "group_id")], stringsAsFactors = FALSE)
    for (m in met_cols) {
      v <- dg[[m]]
      mu  <- stats::median(v, na.rm = TRUE)
      den <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
      if (!is.finite(den) ||
          den == 0)
        den <- stats::IQR(v, na.rm = TRUE) / 1.349
      if (!is.finite(den) ||
          den == 0)
        den <- stats::sd(v, na.rm = TRUE)
      if (!is.finite(den) || den == 0)
        den <- 1
      Zg[[m]] <- (v - mu) / den
    }
    Zg
  }))
  
  # Candidates: flagged non-QC samples only
  z_long <- tidyr::pivot_longer(
    z_tbl,
    cols = all_of(met_cols),
    names_to = "metabolite",
    values_to = "z"
  )
  cand_long <- subset(z_long,
                      group_id != "QC" & is.finite(z) &
                        (z >= z_threshold | z <= -z_threshold))
  
  push_row <- function(rows, r, method, pval, ratio, decision) {
    df <- rows[r, , drop = FALSE]
    df$method <- method
    df$p_value <- pval
    df$test_strength <- ratio
    df$decision <- decision
    df
  }
  
  confirm_rows <- list()
  if (nrow(cand_long) > 0) {
    cand_long <- merge(cand_long, qc_rsd_tbl, by = "metabolite", all.x = TRUE)
    by_key <- interaction(cand_long$group_id, cand_long$metabolite, drop = TRUE)
    split_idx <- split(seq_len(nrow(cand_long)), by_key)
    
    for (idx in split_idx) {
      rows <- cand_long[idx, , drop = FALSE]
      gid  <- rows$group_id[1]
      met  <- rows$metabolite[1]
      
      grp  <- df2[df2$group_id == gid &
                    df2$class != "QC", c("sample", met), drop = FALSE]
      vals <- grp[[met]]
      names(vals) <- grp$sample
      valid <- which(is.finite(vals))
      if (length(valid) < min_group_n) {
        for (r in seq_len(nrow(rows))) {
          confirm_rows[[length(confirm_rows) + 1]] <-  push_row(
            rows,
            r,
            method = NA_character_,
            pval = NA_real_,
            ratio = NA_real_,
            decision = "insufficient_n"
          )
        }
        next
      }
      
      rsd <- rows$qc_rsd[1]
      gate <- ifelse(
        !is.na(rsd) & rsd <= 100 * qc_rsd_stable,
        "stable",
        ifelse(
          !is.na(rsd) &
            rsd > 100 * qc_rsd_unstable,
          "unstable",
          "borderline"
        )
      )
      
      for (r in seq_len(nrow(rows))) {
        s <- rows$sample[r]
        x <- as.numeric(vals[valid])
        names(x) <- names(vals)[valid]
        
        if (!(s %in% names(x))) {
          confirm_rows[[length(confirm_rows) + 1]] <- push_row(
            rows,
            r,
            method = NA_character_,
            pval = NA_real_,
            ratio = NA_real_,
            decision = "no_value"
          )
          next
        }
        if (gate == "unstable") {
          confirm_rows[[length(confirm_rows) + 1]] <- push_row(
            rows,
            r,
            method = NA_character_,
            pval = NA_real_,
            ratio = NA_real_,
            decision = "skip_unstable_qc"
          )
          next
        }
        # stricter z for borderline QC
        z_ok <- abs(rows$z[r]) >= ifelse(gate == "borderline", max(z_threshold, 5), z_threshold)
        if (!z_ok) {
          confirm_rows[[length(confirm_rows) + 1]] <- push_row(
            rows,
            r,
            method = NA_character_,
            pval = NA_real_,
            ratio = NA_real_,
            decision = "below_borderline_z"
          )
          next
        }
        if (!has_out) {
          confirm_rows[[length(confirm_rows) + 1]] <- push_row(
            rows,
            r,
            method = NA_character_,
            pval = NA_real_,
            ratio = NA_real_,
            decision = "outliers_pkg_missing"
          )
          next
        }
        
        # MD contribution as an independent confirm signal
        md_row <- md_tbl[md_tbl$sample == s &
                           md_tbl$group_id == gid, , drop = FALSE]
        md_ok  <- nrow(md_row) > 0 && isTRUE(md_row$flagged[1])
        
        n <- length(x)
        # For small/medium-n tests we need extreme status
        is_min <- which.min(x) == which(names(x) == s)
        is_max <- which.max(x) == which(names(x) == s)
        tied_extreme <- (sum(x == min(x)) > 1) ||
          (sum(x == max(x)) > 1)
        
        
        
        # holders
        meth <- NA_character_
        pval <- NA_real_
        ratio <- NA_real_
        decision <- "retain"
        tt <- NULL
        
        # Large-n: Rosner (ESD)
        if (n > 25) {
          if (has_env) {
            K <- min(max(1, floor(0.1 * n)), 10, n - 2)
            tt <- tryCatch(
              EnvStats::rosnerTest(x, k = K, alpha = alpha),
              error = function(e)
                NULL
            )
            if (is.null(tt)) {
              meth <- "rosner"
              decision <- "test_error"
            } else {
              stats <- tt$all.stats
              obs_col <- grep("^Obs(\\.|_)?Num$", names(stats), value = TRUE)[1]
              r_col   <- grep("^R", names(stats), value = TRUE)[1]
              lam_col <- grep("^lambda",
                              names(stats),
                              ignore.case = TRUE,
                              value = TRUE)[1]
              out_col <- grep("^Outlier$",
                              names(stats),
                              ignore.case = TRUE,
                              value = TRUE)[1]
              pos <- match(s, names(x))
              row_idx <- if (!is.na(obs_col))
                match(pos, stats[[obs_col]])
              else
                NA_integer_
              out_mask <- if (!is.na(out_col) &&
                              is.logical(stats[[out_col]]))
                stats[[out_col]]
              else if (!is.na(out_col))
                stats[[out_col]] %in% c("Yes", "TRUE", "True", "T", "1")
              else
                rep(FALSE, nrow(stats))
              is_out <- is.finite(row_idx) &&
                !is.na(row_idx) && isTRUE(out_mask[row_idx])
              ratio <- if (!is.na(row_idx) &&
                           !is.na(r_col) && !is.na(lam_col)) {
                as.numeric(stats[row_idx, r_col]) / as.numeric(stats[row_idx, lam_col])
              } else
                NA_real_
              meth <- "rosner"
              decision <- if (is_out)
                "confirm"
              else
                "retain"
              pval <- NA_real_
            }
          } else {
            meth <- "rosner"
            decision <- "pkg_missing"
          }
          
        } else {
          # Small/medium-n: Dixon or Grubbs guarded
          use_dixon <- (confirm_method == "dixon") &&
            (n >= 3 && n <= 30)
          if (use_dixon)
            use_dixon <- (is_min || is_max) && !tied_extreme
          
          if (use_dixon) {
            tt <- tryCatch(
              outliers::dixon.test(x, two.sided = TRUE),
              error = function(e)
                NULL
            )
            meth <- if (is.null(tt))
              "grubbs"
            else
              "dixon"
            if (is.null(tt))
              tt <- tryCatch(
                outliers::grubbs.test(x, two.sided = TRUE),
                error = function(e)
                  NULL
              )
          } else {
            if (!(is_min || is_max) || tied_extreme) {
              if (md_ok) {
                # MD override confirms even when Grubbs cannot run
                confirm_rows[[length(confirm_rows) + 1]] <-
                  push_row(
                    rows,
                    r,
                    method = "md_only",
                    pval = NA_real_,
                    ratio = NA_real_,
                    decision = "confirm"
                  )
              } else {
                confirm_rows[[length(confirm_rows) + 1]] <-
                  push_row(
                    rows,
                    r,
                    method = "grubbs",
                    pval = NA_real_,
                    ratio = NA_real_,
                    decision = "not_extreme_for_grubbs"
                  )
              }
              next
            }
            tt <- tryCatch(
              outliers::grubbs.test(x, two.sided = TRUE),
              error = function(e)
                NULL
            )
            meth <- "grubbs"
          }
          
          if (!is.null(tt)) {
            decision <- ifelse(tt$p.value < alpha, "confirm", "retain")
            pval <- as.numeric(tt$p.value)
          } else {
            decision <- "test_error"
          }
        }
        
        # Final decision: test_ok OR md_ok
        if (decision %in% c("confirm", "retain")) {
          test_ok <- (decision == "confirm")
          decision <- if ((test_ok ||
                           md_ok))
            "confirm"
          else
            "retain"
        } else if (md_ok) {
          decision <- "confirm"
          meth <- if (is.na(meth))
            "md_only"
          else
            paste0(meth, "+md")
        }
        
        confirm_rows[[length(confirm_rows) + 1]] <-
          push_row(
            rows,
            r,
            method = meth,
            pval = pval,
            ratio = ratio,
            decision = decision
          )
      }
    }
  }
  
  confirm_tbl <- if (length(confirm_rows))
    do.call(rbind, confirm_rows)
  else
    data.frame(
      sample = character(),
      batch = character(),
      class = character(),
      order = numeric(),
      group_id = character(),
      metabolite = character(),
      z = numeric(),
      qc_rsd = numeric(),
      method = character(),
      p_value = numeric(),
      test_strength = numeric(),
      decision = character(),
      row.names = NULL
    )
  
  list(
    qc_rsd = qc_rsd_tbl,
    sample_md = md_tbl,
    candidate_metabolites = cand_long,
    confirmations = confirm_tbl,
    params = list(
      group_nonqc_by_class = group_nonqc_by_class,
      z_threshold = z_threshold,
      alpha = alpha
    )
  )
}
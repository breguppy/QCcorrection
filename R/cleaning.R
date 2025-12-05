#' Clean data and track replacements
#'
#' @keywords internal
#' @noRd
clean_data <- function(df,
                       sample,
                       batch,
                       class,
                       order,
                       withheld_cols) {
  if (!(batch %in% colnames(df))) {
    df$batch <- "batch1"
  }
  # remove columns that will not be corrected and are not metadata
  df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]
  
  #TODO: keep withheld columns so they may be added back post correction.
  
  # Rename metadata columns
  names(df)[names(df) == sample] <- "sample"
  names(df)[names(df) == batch]  <- "batch"
  names(df)[names(df) == class]  <- "class"
  names(df)[names(df) == order]  <- "order"
  
  # get names of non-numeric columns that are not metadata columns
  non_numeric_cols <- names(df)[sapply(df, function(col) {
    vals <- col[!is.na(col)]
    # TRUE if no value can be converted to numeric
    all(is.na(suppressWarnings(as.numeric(vals))))
  })]
  non_numeric_cols <- setdiff(non_numeric_cols, c("sample", "batch", "class", "order"))
  
  df <- df[, !(names(df) %in% non_numeric_cols)]
  # make sure data is in injection order
  df <- df[order(df$order), ]
  metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
  
  # replace any non-numeric values or exactly 0 values with NA
  # count them as missing values
  repl <- tibble(
    metabolite = metab,
    non_numeric_replaced = 0L,
    zero_replaced = 0L
  )
  for (i in seq_along(metab)) {
    col <- metab[i]
    orig <- df[[col]]
    num  <- suppressWarnings(as.numeric(orig))
    cnt1 <- sum(is.na(num) & !is.na(orig))
    cnt2 <- sum(num == 0, na.rm = TRUE)
    num[num == 0] <- NA
    df[[col]]   <- num
    repl$non_numeric_replaced[i] <- cnt1
    repl$zero_replaced[i]        <- cnt2
  }
  
  # Rename qc samples to be "QC" in class column
  df$class[is.na(df$class)] <- "QC"
  df$class[df$class %in% c("qc", "Qc")] <- "QC"
  
  # make sure data starts and ends with a QC
  if (df$class[1] != "QC") {
    stop("Data sorted by injection order must begin with a QC sample.")
  } else if (df$class[nrow(df)] != "QC") {
    stop("Data sorted by injection order must end with a QC sample.")
  }
  
  # equal columns:
  duplicate_mets <- find_equal_metabolite_cols(df, metab, tolerance = 1e-3)
  correlated_mets <- find_highly_correlated_metabolite_cols(df, metab)
  high_corr_mets <- correlated_mets$high_corr
  print(high_corr_mets)
  
  # Helper to build an unordered pair key
  pair_key <- function(col1, col2) {
    paste(pmin(col1, col2), pmax(col1, col2), sep = "__")
  }
  
  # Only run if both data frames have rows
  if (nrow(duplicate_mets) > 0L && nrow(high_corr_mets) > 0L) {
    dup_keys <- pair_key(duplicate_mets$col1, duplicate_mets$col2)
    cor_keys <- pair_key(high_corr_mets$col1, high_corr_mets$col2)
    
    # Keep only correlated pairs that are not in duplicates
    high_corr_mets <- high_corr_mets[!(cor_keys %in% dup_keys), , drop = FALSE]
  } else {
    high_corr_mets <- high_corr_mets
  }
  
  return(list(
    df = df,
    replacement_counts = repl,
    withheld_cols = withheld_cols,
    non_numeric_cols = non_numeric_cols,
    duplicate_mets = duplicate_mets,
    high_corr_mets = high_corr_mets,
    all_corr_mets = correlated_mets$all_corr
  ))
}

#' Find (nearly) equal columns ignoring NAs
#'
#' @param df   A data frame containing metabolite columns.
#' @param cols Optional character vector of column names to check.
#'             If NULL, all columns in `df` are used.
#' @param ...  Additional arguments passed to `all.equal()`
#'             (e.g. tolerance = 1e-8).
#'
#' @return A data.frame with columns `col1` and `col2` listing ordered pairs
#'         of columns that are equal (per `all.equal`) on all rows where
#'         both columns are non-NA. If no pairs are found, returns an empty
#'         data.frame with the same columns.
#' @keywords internal
#' @noRd
find_equal_metabolite_cols <- function(df, cols = NULL, ...) {
  if (is.null(cols)) {
    cols <- names(df)
  }
  
  cols <- intersect(cols, names(df))
  if (length(cols) < 2L) {
    return(data.frame(col1 = character(0), col2 = character(0)))
  }
  
  # All unique unordered column pairs
  pairs <- utils::combn(cols, 2L, simplify = FALSE)
  
  results <- vector("list", length(pairs))
  keep    <- logical(length(pairs))
  
  for (i in seq_along(pairs)) {
    c1 <- pairs[[i]][1L]
    c2 <- pairs[[i]][2L]
    
    x <- df[[c1]]
    y <- df[[c2]]
    
    # Only compare rows where both are non-NA
    idx <- !is.na(x) & !is.na(y)
    if (!any(idx)) {
      next
    }
    
    # all.equal() returns TRUE or a string; we only want TRUE
    if (isTRUE(all.equal(x[idx], y[idx], ...))) {
      results[[i]] <- data.frame(col1 = c1, col2 = c2, stringsAsFactors = FALSE)
      keep[i] <- TRUE
    }
  }
  
  if (!any(keep)) {
    return(data.frame(col1 = character(0), col2 = character(0)))
  }
  
  do.call(rbind, results[keep])
}

#' Find highly correlated metabolite columns ignoring NAs
#'
#' @param df            Data frame containing metabolite columns.
#' @param cols          Optional vector of column names to check.
#' @param method        Correlation method passed to stats::cor().
#' @param cor_threshold Minimum absolute correlation to flag as "highly correlated".
#' @param min_complete  Minimum paired non-NA observations required.
#'
#' @return A list with:
#'   - high_corr : data.frame of pairs passing the threshold
#'   - all_corr  : data.frame of all computed correlations, sorted descending
#'
#' @keywords internal
#' @noRd
find_highly_correlated_metabolite_cols <- function(df,
                                                   cols = NULL,
                                                   method = "pearson",
                                                   cor_threshold = 0.995,
                                                   min_complete = 3L) {
  
  # Determine candidate columns
  if (is.null(cols)) cols <- names(df)
  cols <- intersect(cols, names(df))
  
  # Restrict to numeric columns
  is_num <- vapply(df[cols], is.numeric, logical(1L))
  cols   <- cols[is_num]
  
  if (length(cols) < 2L) {
    empty <- data.frame(col1=character(0), col2=character(0),
                        cor=numeric(0), n_complete=integer(0))
    return(list(high_corr = empty, all_corr = empty))
  }
  
  # All unordered pairs
  pairs <- utils::combn(cols, 2L, simplify = FALSE)
  
  high_results <- vector("list", length(pairs))
  all_results  <- vector("list", length(pairs))
  
  keep_high <- logical(length(pairs))
  keep_all  <- logical(length(pairs))
  
  for (i in seq_along(pairs)) {
    c1 <- pairs[[i]][1L]
    c2 <- pairs[[i]][2L]
    
    x <- df[[c1]]
    y <- df[[c2]]
    
    idx <- is.finite(x) & is.finite(y)
    n_complete <- sum(idx)
    
    if (n_complete < min_complete) next
    
    # Compute correlation
    cor_val <- suppressWarnings(
      stats::cor(x[idx], y[idx], method = method)
    )
    
    if (is.na(cor_val)) next
    
    # Store in 'all correlations'
    all_results[[i]] <- data.frame(
      col1 = c1, col2 = c2,
      cor = unname(cor_val),
      n_complete = n_complete,
      stringsAsFactors = FALSE
    )
    keep_all[i] <- TRUE
    
    # Store only if above threshold
    if (abs(cor_val) >= cor_threshold) {
      high_results[[i]] <- all_results[[i]]
      keep_high[i] <- TRUE
    }
  }
  
  # Construct outputs
  all_corr <- if (any(keep_all)) {
    out <- do.call(rbind, all_results[keep_all])
    out[order(-abs(out$cor)), ]   # Sort descending by absolute correlation
  } else {
    data.frame(col1=character(0), col2=character(0),
               cor=numeric(0), n_complete=integer(0))
  }
  
  high_corr <- if (any(keep_high)) {
    do.call(rbind, high_results[keep_high])
  } else {
    data.frame(col1=character(0), col2=character(0),
               cor=numeric(0), n_complete=integer(0))
  }
  
  list(
    high_corr = high_corr,
    all_corr  = all_corr
  )
}


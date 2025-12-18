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
  
  return(list(
    df = df,
    replacement_counts = repl,
    withheld_cols = withheld_cols,
    non_numeric_cols = non_numeric_cols,
    duplicate_mets = duplicate_mets
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
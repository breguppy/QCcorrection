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
  # remove columns that will not be corrected and are not metadata
  df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]
  
  #TODO: keep withheld columns so they may be added back post correction.
  
  # Rename metadata columns
  names(df)[names(df) == sample] <- "sample"
  names(df)[names(df) == batch]  <- "batch"
  names(df)[names(df) == class]  <- "class"
  names(df)[names(df) == order]  <- "order"
  
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
  
  return(list(
    df = df,
    replacement_counts = repl,
    withheld_cols = withheld_cols
  ))
}
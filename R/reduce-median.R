#' Median across model outputs with key alignment
#' @keywords internal
#' @noRd
.median_across_models <- function(df_list,
                                  metadata_cols = c("sample","batch","class","order")) {
  stopifnot(length(df_list) >= 1L)
  
  # schema check
  ref_cols <- names(df_list[[1]])
  for (i in seq_along(df_list)) {
    if (!setequal(names(df_list[[i]]), ref_cols))
      stop("Model data frames have different columns.")
  }
  
  metab_cols <- setdiff(ref_cols, metadata_cols)
  if (!length(metab_cols)) stop("No metabolite columns detected.")
  
  # enforce numeric metabolite columns
  for (i in seq_along(df_list)) {
    nn <- !vapply(df_list[[i]][metab_cols], is.numeric, TRUE)
    if (any(nn)) df_list[[i]][metab_cols[nn]] <- lapply(df_list[[i]][metab_cols[nn]], as.numeric)
  }
  
  make_key <- function(d) {
    if (any(!metadata_cols %in% names(d)))
      stop("Missing metadata columns in a model output.")
    if (any(duplicated(d[metadata_cols])))
      stop("Duplicate key rows found in a model output.")
    do.call(paste, c(d[metadata_cols], sep = "\r"))
  }
  
  keys_list   <- lapply(df_list, make_key)
  common_keys <- Reduce(intersect, keys_list)
  if (!length(common_keys)) stop("No overlapping rows across model outputs.")
  
  # align each df to common keys
  align_to <- function(d, keys) {
    k <- make_key(d)
    idx <- match(keys, k)
    if (anyNA(idx)) stop("Internal alignment failure.")
    d[idx, c(metadata_cols, metab_cols), drop = FALSE]
  }
  aligned <- lapply(df_list, align_to, keys = common_keys)
  
  # stack metabolite matrices: row x metab x models
  mats <- lapply(aligned, function(d) as.matrix(d[metab_cols]))
  arr  <- simplify2array(mats)
  med  <- apply(arr, c(1, 2), stats::median, na.rm = TRUE)
  
  out <- cbind(aligned[[1]][metadata_cols], as.data.frame(med))
  names(out) <- c(metadata_cols, metab_cols)
  
  # stable order
  if (all(c("batch","order") %in% names(out))) {
    out <- out[order(out$batch, out$order), , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}
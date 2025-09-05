#' RSD functions
#'
#' @keywords internal
#' @noRd
metabolite_rsd <- function(df,
                           metadata_cols = c("sample", "batch", "class", "order")) {
  nm <- names(df)
  md_idx <- tolower(nm) %in% tolower(metadata_cols)
  class_col <- nm[tolower(nm) == "class"]
  if (!length(class_col))
    stop("Expected a 'class' column (any case).")
  class_col <- class_col[1]
  metab_cols <- nm[!md_idx]
  is_num <- vapply(df[metab_cols], is.numeric, logical(1))
  metab_cols <- metab_cols[is_num]
  if (!length(metab_cols))
    stop("No numeric metabolite columns detected.")
  
  # Separate QC and non-QC samples
  qc_df <- df[df[[class_col]] == "QC", metab_cols, drop = FALSE]
  nonqc_df <- df[df[[class_col]] != "QC", metab_cols, drop = FALSE]
  
  # RSD function ignoring NA
  rsd_fun <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- stats::sd(x, na.rm = TRUE)
    if (mu == 0 || !is.finite(mu))
      return(NA_real_)
    return(100 * sigma / mu)
  }
  
  data.frame(
    Metabolite = metab_cols,
    RSD_QC     = vapply(qc_df, rsd_fun, numeric(1)),
    RSD_NonQC  = vapply(nonqc_df, rsd_fun, numeric(1)),
    check.names = FALSE
  )
}


class_metabolite_rsd <- function(df,
                                 metadata_cols = c("sample", "batch", "class", "order")) {
  # Get metabolite columns
  metab_cols <- setdiff(names(df), metadata_cols)
  
  long_df <- df %>%
    pivot_longer(
      cols = all_of(metab_cols),
      names_to = "Metabolite",
      values_to = "Value"
    )
  
  # Group by class and metabolite and compute RSD
  rsd_df <- long_df %>%
    group_by(class, Metabolite) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      RSD = ifelse(Mean == 0, NA, SD / Mean * 100),
      .groups = "drop"
    )
  
  return(rsd_df)
}
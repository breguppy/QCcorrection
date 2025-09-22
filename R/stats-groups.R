#' Group stats functions
#'
#' @keywords internal
#' @noRd
group_stats <- function(df) {
  # Get metabolite columns, rename class to Group, remove batch and order
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  names(df)[names(df) == "class"] <- "Group"
  df$batch <- NULL
  df$order <- NULL
  
  # Compute statistics for each group
  group_dfs = list()
  group_stats_dfs = list()
  group_names <- unique(df$Group)
  for (group_name in group_names) {
    group_df <- df %>% filter(Group == group_name)
    group_dfs[[group_name]] <- group_df
    
    means <- summarise(group_df, across(all_of(metab_cols), ~ mean(., na.rm = TRUE))) %>%
      mutate(` ` = "Mean")
    ses <- summarise(group_df, across(all_of(metab_cols), ~ sd(., na.rm = TRUE) / sqrt(sum(!is.na(
      .
    ))))) %>%
      mutate(` ` = "SE")
    cvs <- summarise(group_df, across(
      all_of(metab_cols),
      ~ sd(., na.rm = TRUE) / mean(., na.rm = TRUE)
    )) %>%
      mutate(` ` = "CV")
    
    # Bind summary rows
    group_stats_df <- bind_rows(means, ses, cvs) %>%
      select(` `, all_of(metab_cols))
    group_stats_dfs[[group_name]] <- group_stats_df
  }
  
  return(list(group_dfs = group_dfs, group_stats_dfs = group_stats_dfs))
}

fold_changes <- function(df, control_mean) {
  # Get metabolite columns
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  
  # Check formatting and ensure there is a control mean
  if (nrow(control_mean) != 1) {
    stop("control_mean must contain exactly one row with metabolite means.")
  }
  
  fold_change <- df
  
  # Divide each metabolite column by the corresponding control mean
  for (col in metab_cols) {
    fold_change[[col]] <- fold_change[[col]] / control_mean[[col]][1]
  }
  
  return(fold_change)
}
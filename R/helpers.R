# Helper function for app

library(statTarget)

# Remove metabolites based on Frule
filter_data <- function(df, metabolite_cols, Frule){
  # remove any metabolite column that has more than Frule% missing values
  
  # Compute percentage of missing values per metabolite column
  missing_pct <- sapply(df[metabolite_cols], function(col) {
    mean(is.na(col)) * 100
  })
  
  # Keep only columns with missing percentage <= Frule
  keep_cols <- metabolite_cols[missing_pct <= Frule]
  removed_cols <- setdiff(metabolite_cols, keep_cols)
  
  # Return df with only the retained metabolite columns
  df_filtered <- df[, c(setdiff(names(df), metabolite_cols), keep_cols)]
  
  return(list(
    df_filtered = df_filtered,
    removed_cols = removed_cols
    ))
}

# Impute missing values based on imputeM
impute_missing <- function(df, metabolite_cols, imputeM) {
  # impute missing values based on the imputation method
  
}

# prep data for correction


# Correct data
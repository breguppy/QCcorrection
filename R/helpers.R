# Helper function for app

library(statTarget)
library(dplyr)

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
impute_missing <- function(df, metabolite_cols, imputeM, sample_class) {
  
  imputed_df <- df
  
  # impute missing values with column median
  if (imputeM == "median"){
    imputed_df[metabolite_cols] <- lapply(imputed_df[metabolite_cols], function(col) {
      col[is.na(col)] <- median(col, na.rm = TRUE)
      col
    })
  }
  # impute missing values with column mean
  else if (imputeM == "mean"){
    imputed_df[metabolite_cols] <- lapply(imputed_df[metabolite_cols], function(col) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      col
    })
  }
  # Impute with sample_class median for each column
  else if (imputeM == "class_median"){
    imputed_df <- imputed_df %>%
      group_by(.data[[sample_class]]) %>%
      mutate(across(all_of(metabolite_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with sample_class mean for each column
  else if (imputeM == "class_mean"){
    imputed_df <- imputed_df %>%
      group_by(.data[[sample_class]]) %>%
      mutate(across(all_of(metabolite_cols), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with KNN
  else if (imputeM == "KNN"){
    metab_matrix <- as.matrix(imputed_df[metabolite_cols])
    transposed <- t(metab_matrix)
    
    knn_result <- impute::impute.knn(transposed, rowmax = 0.99, colmax = 0.99, 
                                     maxp = 15000)
    
    imputed_matrix <- t(knn_result$data)
    imputed_df[metabolite_cols] <- as.data.frame(imputed_matrix)
  }
  # impute with minimum column value
  else if (imputeM == "min"){
    imputed_df[metabolite_cols] <- lapply(imputed_df[metabolite_cols], function(col) {
      col[is.na(col)] <- min(col, na.rm = TRUE)
      col
    })
  }
  # impute with half minimum column value
  else if (imputeM == "minHalf"){
    imputed_df[metabolite_cols] <- lapply(imputed_df[metabolite_cols], function(col) {
      col[is.na(col)] <- 0.5 * min(col, na.rm = TRUE)
      col
    })
  }
  # impute with 0
  else if (imputeM == "zero"){
    imputed_df[metabolite_cols] <- lapply(imputed_df[metabolite_cols], function(col) {
      col[is.na(col)] <- 0
      col
    })
  }
  
  return(imputed_df)
}

# prep data for correction


# Correct data
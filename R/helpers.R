# Helper function for app

library(randomForest)
library(dplyr)

# Remove metabolites based on Frule
filter_data <- function(df, metab_cols, Frule){
  # remove any metabolite column that has more than Frule% missing values
  
  # Compute percentage of missing values per metabolite column
  missing_pct <- sapply(df[metab_cols], function(col) {
    mean(is.na(col)) * 100
  })
  
  # Keep only columns with missing percentage <= Frule
  keep_cols <- metab_cols[missing_pct <= Frule]
  removed_cols <- setdiff(metab_cols, keep_cols)
  
  # Return df with only the retained metabolite columns
  df_filtered <- df[, c(setdiff(names(df), metab_cols), keep_cols)]
  
  return(list(
    df_filtered = df_filtered,
    removed_cols = removed_cols
    ))
}

# Impute missing values based on imputeM
impute_missing <- function(df, metab_cols, imputeM, sample_class) {
  n_missv <- sum(is.na(df[, metab_cols]))
  imputed_df <- df
  impute_str <- imputeM
  # impute missing values with column median
  if (imputeM == "median"){
    impute_str <- "metabolite median"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- median(col, na.rm = TRUE)
      col
    })
  }
  # impute missing values with column mean
  else if (imputeM == "mean"){
    impute_str <- "metabolite mean"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      col
    })
  }
  # Impute with sample_class median for each column
  else if (imputeM == "class_median"){
    impute_str <- "class-metabolite median"
    imputed_df <- imputed_df %>%
      group_by(.data[[sample_class]]) %>%
      mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with sample_class mean for each column
  else if (imputeM == "class_mean"){
    impute_str <- "class-metabolite median"
    imputed_df <- imputed_df %>%
      group_by(.data[[sample_class]]) %>%
      mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with KNN
  else if (imputeM == "KNN"){
    metab_matrix <- as.matrix(imputed_df[metab_cols])
    transposed <- t(metab_matrix)
    knn_result <- impute::impute.knn(transposed, rowmax = 0.99, colmax = 0.99, 
                                     maxp = 15000)
    imputed_matrix <- t(knn_result$data)
    imputed_df[metab_cols] <- as.data.frame(imputed_matrix)
  }
  # impute with minimum column value
  else if (imputeM == "min"){
    impute_str <- "metabolite minimum"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- min(col, na.rm = TRUE)
      col
    })
  }
  # impute with half minimum column value
  else if (imputeM == "minHalf"){
    impute_str <- "half the metabolite minimum"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- 0.5 * min(col, na.rm = TRUE)
      col
    })
  }
  # impute with 0
  else if (imputeM == "zero"){
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- 0
      col
    })
  }
  
  return(list(
    df_imputed = imputed_df,
    n_missv = n_missv,
    impute_str = impute_str
  ))
}

# QCRFSC (statTarget style)
qc_rfsc_correction <- function(df, metab_cols, ntree = 500, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  df_corrected <- df
  
  pb <- txtProgressBar(min = 0, max = length(metab_cols), style = 3)
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # QC samples: class is NA
    qc_idx <- which(is.na(df$class))
    non_qc_idx <- which(!is.na(df$class))
    
    if (length(qc_idx) < 5) {
      warning(paste("Skipping", metab, "- too few QC samples."))
      next
    }
    
    # Build random forest using QC samples
    rf_model <- randomForest(
      x = data.frame(order = df$order[qc_idx]),
      y = df[[metab]][qc_idx],
      ntree = ntree
    )
    
    # Predict values for all samples using their order
    predicted <- predict(rf_model, newdata = data.frame(order = df$order))
    
    # Apply correction
    df_corrected[[metab]] <- df[[metab]] / predicted
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(df_corrected)
}

# Correct data
correct_data <- function(df, metab_cols, MLmethod){
  
  if (MLmethod == "QCRFSC"){
    cor_str <- "QC Random Forest Signal Correction (statTarget style)"
    # Get QC indexes from class column.
    df_corrected <- qc_rfsc_correction(df, metab_cols, ntree = 500, seed = 123)
    # Build random forest for each metabolite based on QCs
    # Apply correction
  }
  
  return(list(
    df_corrected = df_corrected,
    cor_str = cor_str
  ))
}
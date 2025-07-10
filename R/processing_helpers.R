# Processing helpers

library(randomForest)
library(dplyr)
library(purrr)
library(stats)

#–– Data cleaning helpers ––#

# Clean & track replacements
cleanData <- function(df, sample, batch, class, order) {
  names(df)[names(df) == sample] <- "sample"
  names(df)[names(df) == batch]  <- "batch"
  names(df)[names(df) == class]  <- "class"
  names(df)[names(df) == order]  <- "order"
  
  df <- df[order(df$order), ]
  metab <- setdiff(names(df), c("sample", "batch", "class", "order"))
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
  
  df$class[is.na(df$class)] <- "QC"
  df$class[df$class %in% c("qc", "Qc")] <- "QC"
  
  if (df$class[1] != "QC") {
    stop("Data sorted by injection order must begin  with a QC sample.")
  } else if (df$class[nrow(df)] != "QC") {
    stop("Data sorted by injection order must end with a QC sample.")
  }
  
  list(df = df, replacement_counts = repl)
}

# Remove metabolites based on Frule
filter_data <- function(df, metab_cols, Frule) {
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
  
  return(list(df_filtered = df_filtered, removed_cols = removed_cols))
}


# Impute missing values based on imputeM
impute_missing <- function(df, metab_cols, imputeM) {
  n_missv <- sum(is.na(df[, metab_cols]))
  imputed_df <- df
  impute_str <- imputeM
  # impute missing values with column median
  if (imputeM == "median") {
    impute_str <- "metabolite median"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- median(col, na.rm = TRUE)
      col
    })
  }
  # impute missing values with column mean
  else if (imputeM == "mean") {
    impute_str <- "metabolite mean"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      col
    })
  }
  # Impute with sample_class median for each column
  else if (imputeM == "class_median") {
    impute_str <- "class-metabolite median"
    imputed_df <- imputed_df %>%
      group_by(.data[["class"]]) %>%
      mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with sample_class mean for each column
  else if (imputeM == "class_mean") {
    impute_str <- "class-metabolite median"
    imputed_df <- imputed_df %>%
      group_by(.data[["class"]]) %>%
      mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
      ungroup()
  }
  # impute with KNN
  else if (imputeM == "KNN") {
    metab_matrix <- as.matrix(imputed_df[metab_cols])
    transposed <- t(metab_matrix)
    knn_result <- impute::impute.knn(transposed,
                                     rowmax = 0.99,
                                     colmax = 0.99,
                                     maxp = 15000)
    imputed_matrix <- t(knn_result$data)
    imputed_df[metab_cols] <- as.data.frame(imputed_matrix)
  }
  # impute with minimum column value
  else if (imputeM == "min") {
    impute_str <- "metabolite minimum"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- min(col, na.rm = TRUE)
      col
    })
  }
  # impute with half minimum column value
  else if (imputeM == "minHalf") {
    impute_str <- "half the metabolite minimum"
    imputed_df[metab_cols] <- lapply(imputed_df[metab_cols], function(col) {
      col[is.na(col)] <- 0.5 * min(col, na.rm = TRUE)
      col
    })
  }
  # impute with 0
  else if (imputeM == "zero") {
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


#–– Data correction helpers ––#

# 3 by 500 random forest correction
rf_correction <- function(df,
                          metab_cols,
                          ntree = 500,
                          seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  df_corrected <- df
  
  pb <- txtProgressBar(min = 0,
                       max = length(metab_cols),
                       style = 3)
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # QC samples: class is NA
    qc_idx <- which(df$class == "QC")
    non_qc_idx <- which(df$class != "QC")
    
    if (length(qc_idx) < 5) {
      warning(paste("Skipping", metab, "- too few QC samples."))
      next
    }
    
    # Build random forest using QC samples
    rf_model <- randomForest::randomForest(x = data.frame(order = df$order[qc_idx]),
                                           y = df[[metab]][qc_idx],
                                           ntree = ntree)
    
    # Predict values for all samples using their order
    predicted <- predict(rf_model, newdata = data.frame(order = df$order))
    
    # Apply correction
    df_corrected[[metab]] <- df[[metab]] / predicted
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(df_corrected)
}

compute_median_dataframe <- function(df_list,
                                     metadata_cols = c("Sample", "Batch", "Order", "Treatment")) {
  # Get metabolite column names from the first dataframe
  all_cols <- colnames(df_list[[1]])
  metabolite_cols <- setdiff(all_cols, metadata_cols)
  
  # Add a composite key for alignment
  df_list_keyed <- lapply(df_list, function(df) {
    df$key <- do.call(paste, c(df[metadata_cols], sep = "||"))
    df
  })
  
  # Sort by the composite key
  df_list_sorted <- lapply(df_list_keyed, function(df)
    df[order(df$key), ])
  
  # Sanity check: all keys must align
  key_ref <- df_list_sorted[[1]]$key
  if (!all(sapply(df_list_sorted, function(df)
    all(df$key == key_ref)))) {
    stop("Composite keys are not aligned across all data frames.")
  }
  
  # Extract aligned metabolite matrices
  metabolite_array <- array(sapply(df_list_sorted, function(df)
    as.matrix(df[, metabolite_cols])),
    dim = c(
      nrow(df_list_sorted[[1]]),
      length(metabolite_cols),
      length(df_list_sorted)
    ))
  
  # Compute median along third dimension
  median_matrix <- apply(metabolite_array, c(1, 2), median, na.rm = TRUE)
  
  # Reconstruct result
  median_df <- df_list_sorted[[1]][, metadata_cols]
  median_df <- cbind(median_df, as.data.frame(median_matrix))
  colnames(median_df)[-(1:length(metadata_cols))] <- metabolite_cols
  
  median_df <- median_df[order(median_df$order), ]
  
  return(median_df)
}

loess_correction <- function(df,
                             metab_cols,
                             degree = 2,
                             span = 0.75) {
  df_corrected <- df
  
  pb <- txtProgressBar(min = 0,
                       max = length(metab_cols),
                       style = 3)
  
  for (i in seq_along(metab_cols)) {
    metab <- metab_cols[i]
    
    # QC samples
    qc_idx <- which(df$class == "QC")
    
    dat <- data.frame(intensity = df[[metab]][qc_idx], order     = df$order[qc_idx])
    
    # Fit loess (intensity ~ order)
    loess_model <- stats::loess(
      intensity ~ order,
      data   = dat,
      span   = span,
      degree = degree
    )
    
    
    # Predict values for all samples using their order
    predicted <- stats::predict(loess_model, newdata = data.frame(order = df$order))
    
    # Apply correction
    df_corrected[[metab]] <- df[[metab]] / predicted
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(df_corrected)
}
bw_loess_correction <- function(df,
                                metab_cols,
                                degree = 2,
                                span = 0.75) {
  df_corrected <- df
  
  batches <- unique(df$batch)
  
  pb <- txtProgressBar(
    min = 0,
    max = length(metab_cols) * length(batches),
    style = 3
  )
  progress_counter <- 0
  for (metab  in metab_cols) {
    for (b in batches) {
      batch_df <- df[df$batch == b, ]
      batch_idx <- which(df$batch == b)
      
      # QC samples in this batch
      qc_idx <- which(batch_df$class == "QC")
      
      if (length(qc_idx) < 5) {
        warning(
          sprintf(
            "Skipping batch '%s' for metabolite '%s': not enough QC samples",
            b,
            metab
          )
        )
        progress_counter <- progress_counter + 1
        setTxtProgressBar(pb, progress_counter)
        next
      }
      
      dat <- data.frame(intensity = batch_df[[metab]][qc_idx], order     = batch_df$order[qc_idx])
      
      # Fit loess (intensity ~ order)
      loess_model <- stats::loess(
        intensity ~ order,
        data   = dat,
        span   = span,
        degree = degree
      )
      
      # Predict on all samples in this batch
      predicted <- stats::predict(loess_model, newdata = data.frame(order = batch_df$order))
      
      # Apply correction only to this batch
      df_corrected[[metab]][batch_idx] <- df[[metab]][batch_idx] / predicted
      
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
  }
  close(pb)
  return(df_corrected)
}

# Correct data
correct_data <- function(df, metab_cols, corMethod) {
  if (corMethod == "RF") {
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (rf_correction(df, metab_cols, ntree = 500, seed = seed))
      
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- compute_median_dataframe(df_list, metadata_cols)
  } else if (corMethod == "LOESS") {
    df_corrected <- loess_correction(df, metab_cols)
  } else if (corMethod == "BW-RF") {
    seeds <- c(42, 31416, 272)
    #TODO: <- bw_rf_correction(df, metab_cols, ntree = 500, seed = seed)
  } else if (corMethod == "BW_LOESS") {
    df_corrected <- bw_loess_correction(df, metab_cols)
  }
  
  return(df_corrected)
}

remove_imputed_from_corrected <- function(raw_df, corrected_df) {
  # Ensure both data frames are the same shape
  if (!all(dim(raw_df) == dim(corrected_df))) {
    stop("Both data frames must have the same dimensions.")
  }
  
  # Return a new corrected_df with values removed where raw_df is NA
  corrected_df[is.na(raw_df)] <- NA
  return(corrected_df)
}

# compute metabolite RSD
metabolite_rsd <- function(df,
                           metadata_cols = c("sample", "batch", "class", "order")) {
  metab_cols <- setdiff(names(df), metadata_cols)
  
  # Separate QC and non-QC samples
  qc_df <- df[df$class == "QC", metab_cols, drop = FALSE]
  nonqc_df <- df[df$class != "QC", metab_cols, drop = FALSE]
  
  # RSD function ignoring NA
  rsd_fun <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    if (mu == 0 || is.na(mu))
      return(NA)
    return(100 * sigma / mu)
  }
  
  # Compute RSDs
  qc_rsd <- sapply(qc_df, rsd_fun)
  nonqc_rsd <- sapply(nonqc_df, rsd_fun)
  
  # Return as data frame
  result <- data.frame(
    Metabolite = names(qc_rsd),
    RSD_QC = qc_rsd,
    RSD_NonQC = nonqc_rsd,
    row.names = NULL,
    check.names = FALSE
  )
  
  return(result)
}

# compute metabolite RSD
class_metabolite_rsd <- function(df,
                                 metadata_cols = c("sample", "batch", "class", "order")) {
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

# Filter metabolites based on rsd_cutoff
rsd_filter <- function(df,
                       rsd_cutoff,
                       metadata_cols = c("sample", "batch", "class", "order")) {
  # Compute RSD
  rsd_df <- metabolite_rsd(df, metadata_cols)
  
  # Identify which metabolites to keep and remove
  keep_metabolites <- rsd_df$Metabolite[is.na(rsd_df$RSD_QC) |
                                          rsd_df$RSD_QC <= rsd_cutoff]
  remove_metabolites <- rsd_df$Metabolite[!is.na(rsd_df$RSD_QC) &
                                            rsd_df$RSD_QC > rsd_cutoff]
  
  # Columns to retain in filtered data
  final_cols <- c(metadata_cols, keep_metabolites)
  filtered_df = df[, final_cols, drop = FALSE]
  
  # Return a list with the filtered data and removed metabolites
  return(list(filtered_df = filtered_df, removed_metabolites = remove_metabolites))
}

total_ratio_normalization <- function(df, metab_cols) {
  metab_data <- df[, metab_cols, drop = FALSE]
  
  # sum metab_cols values in each row (sample).
  row_sums <- rowSums(metab_data, na.rm = TRUE)
  # compute number of non-missing values in each row (sample).
  non_missing_counts <- rowSums(!is.na(metab_data))
  
  # determine row ratio = (number of non-missing) / (row sum)
  ratios <- non_missing_counts / row_sums
  
  # multiply each row by ratio.
  trn_data <- sweep(metab_data, 1, ratios, FUN = "*")
  
  df[, metab_cols] <- trn_data
  
  return(df)
}


transform_data <- function(df, transform, withheld_cols) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  
  transformed_df <- df
  
  if (transform == "log2") {
    transformed_df[metab_cols] <- log(transformed_df[metab_cols], base = 2)
  } else if (transform == "TRN") {
    metab_cols <- setdiff(metab_cols, withheld_cols)
    transformed_df <- total_ratio_normalization(transformed_df, metab_cols)
  }
  return (transformed_df)
}
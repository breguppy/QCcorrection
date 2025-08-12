# Processing helpers

library(randomForest)
library(dplyr)
library(purrr)
library(stats)

#–– Data cleaning helpers ––#
read_raw_data <- function(file_path) {
  file_ext <- tools::file_ext(file_path)
  
  df <- switch(
    tolower(file_ext),
    "csv" = read.csv(file_path, header = TRUE, check.names = FALSE),
    "xls" = read_excel(file_path),
    "xlsx" = read_excel(file_path),
    stop(
      "Unsupported file type. Please upload a .csv, .xls, or .xlsx file."
    )
  )
  
  return(df)
}

# Clean & track replacements
cleanData <- function(df,
                      sample,
                      batch,
                      class,
                      order,
                      withheld_cols) {
  # Remove columns that will not be corrected and are not metadata
  df <- df[, setdiff(names(df), withheld_cols), drop = FALSE]
  
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

#-- Remove metabolites based on Frule
filter_data <- function(df, metab_cols, Frule) {
  # get metadata columns
  meta_cols <- setdiff(names(df), metab_cols)
  
  # Compute percentage of missing values per metabolite column
  missing_pct <- sapply(df[metab_cols], function(col) {
    mean(is.na(col)) * 100
  })
  
  # Keep only columns with missing percentage <= Frule
  mv_keep_cols <- metab_cols[missing_pct <= Frule]
  # list columns removed due to missing value %
  mv_removed_cols <- setdiff(metab_cols, mv_keep_cols)
  
  # filter data by metabolite missing value
  df_filtered <- df[, c(meta_cols, mv_keep_cols)]
  
  return(list(
    df = df_filtered,
    Frule = Frule,
    mv_removed_cols = mv_removed_cols
  ))
}


#-- Impute missing values based on user selected imputation method
impute_missing <- function(df, metab_cols, qcImputeM, samImputeM) {
  n_missv <- sum(is.na(df[, metab_cols]))
  imputed_df <- df
  
  # helper function to apply an imputation strategy to a subset of data
  apply_impute <- function(sub_df, method) {
    if (method == "nothing_to_impute") {
      sub_df <- sub_df
      str <- "Nothing to Impute"
    } else if (method == "median") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- median(col, na.rm = TRUE)
        col
      })
      str <- "Metabolite Median"
    } else if (method == "mean") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
      str <- "Metabolite Mean"
    } else if (method == "class_median") {
      sub_df <- sub_df %>%
        group_by(.data[["class"]]) %>%
        mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
        ungroup()
      str <- "Class-Metabolite Median"
    } else if (method == "class_mean") {
      sub_df <- sub_df %>%
        group_by(.data[["class"]]) %>%
        mutate(across(all_of(metab_cols), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
        ungroup()
      str <- "Class-Metabolite Mean"
    } else if (method == "min") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- min(col, na.rm = TRUE)
        col
      })
      str <- "Minimum Metabolite Value"
    } else if (method == "minHalf") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- 0.5 * min(col, na.rm = TRUE)
        col
      })
      str <- "Half the Minimum Metabolite Value"
    } else if (method == "zero") {
      sub_df[metab_cols] <- lapply(sub_df[metab_cols], function(col) {
        col[is.na(col)] <- 0
        col
      })
      str <- "Zero"
    } else if (method == "KNN") {
      metab_matrix <- as.matrix(sub_df[metab_cols])
      transposed <- t(metab_matrix)
      knn_result <- impute::impute.knn(transposed,
                                       rowmax = 0.99,
                                       colmax = 0.99,
                                       maxp = 15000)
      imputed_matrix <- t(knn_result$data)
      sub_df[metab_cols] <- as.data.frame(imputed_matrix)
      str <- "KNN"
    }
    return(list(sub_df = sub_df, str = str))
  }
  
  # Identify QC and Sample rows
  is_qc <- df$class == "QC"
  qc_df <- imputed_df[is_qc, ]
  sam_df <- imputed_df[!is_qc, ]
  
  # Apply user chosen imputation methods
  qc_imputed <- apply_impute(qc_df, qcImputeM)
  qc_df <- qc_imputed$sub_df
  qc_str <- qc_imputed$str
  sam_imputed <- apply_impute(sam_df, samImputeM)
  sam_df <- sam_imputed$sub_df
  sam_str <- sam_imputed$str
  
  # Combine the results and make sure its in order
  imputed_df <- bind_rows(qc_df, sam_df) %>%
    arrange(.data[["order"]])
  
  return(list(
    df = imputed_df,
    qc_str = qc_str,
    sam_str = sam_str,
    n_missv = n_missv
  ))
}


#–– Data correction helpers ––#

#-- 3 by 500 random forest correction
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
    
    # increase printed progress bar
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

bw_rf_correction <- function(df,
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
    
    for (batch_id in unique(df$batch)) {
      batch_df <- df[df$batch == batch_id, ]
      batch_idx <- which(df$batch == batch_id)
      qc_idx <- which(batch_df$class == "QC")
      
      # Skip if not enough QCs in this batch
      if (length(qc_idx) < 5) {
        warning(sprintf(
          "Skipping %s in batch %s - too few QC samples.",
          metab,
          batch_id
        ))
        next
      }
      
      # Fit random forest using QC samples in this batch
      rf_model <- randomForest::randomForest(
        x = data.frame(order = batch_df$order[qc_idx]),
        y = batch_df[[metab]][qc_idx],
        ntree = ntree
      )
      
      # Predict all values in the batch
      predicted <- predict(rf_model, newdata = data.frame(order = batch_df$order))
      
      # Apply correction
      df_corrected[[metab]][batch_idx] <- df[[metab]][batch_idx] / predicted
    }
    
    # increase printed progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(df_corrected)
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
    
    # increase printed progress bar
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
      
      # increase printed progress bar
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
    correction_str <- "Random Forest"
    parameters <- "3 models with seeds 42, 31416, 272. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (rf_correction(df, metab_cols, ntree = 500, seed = seed))
      
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- compute_median_dataframe(df_list, metadata_cols)
  } else if (corMethod == "LOESS") {
    correction_str <- "LOESS"
    parameters <- "degree = 2 and span = 0.75."
    df_corrected <- loess_correction(df, metab_cols)
  } else if (corMethod == "BW_RF") {
    correction_str <- "Batchwise Random Forest"
    parameters <- "3 models with seeds 42, 31416, 272. Final corrected data is the median value of the 3 models."
    seeds <- c(42, 31416, 272)
    df_list <- lapply(seeds, function(seed) {
      return (bw_rf_correction(df, metab_cols, ntree = 500, seed = seed))
    })
    metadata_cols <- setdiff(colnames(df), metab_cols)
    df_corrected <- compute_median_dataframe(df_list, metadata_cols)
  } else if (corMethod == "BW_LOESS") {
    correction_str <- "Batchwise LOESS"
    parameters <- "degree = 2 and span = 0.75."
    df_corrected <- bw_loess_correction(df, metab_cols)
  }
  
  return(list(
    df = df_corrected,
    str = correction_str,
    parameters =  parameters
  ))
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
  nm <- names(df)
  md_idx <- tolower(nm) %in% tolower(metadata_cols)
  class_col <- nm[tolower(nm) == "class"]
  if (!length(class_col)) stop("Expected a 'class' column (any case).")
  class_col <- class_col[1]
  metab_cols <- nm[!md_idx]
  is_num <- vapply(df[metab_cols], is.numeric, logical(1))
  metab_cols <- metab_cols[is_num]
  if (!length(metab_cols)) stop("No numeric metabolite columns detected.")
  
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
  return(
    list(
      df = filtered_df,
      rsd_cutoff = rsd_cutoff,
      removed_metabolites = remove_metabolites
    )
  )
}

total_ratio_norm <- function(df, metab_cols) {
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


transform_data <- function(df, transform, withheld_cols, ex_ISTD = TRUE) {
  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols <- setdiff(names(df), meta_cols)
  if (ex_ISTD) {
    # Get column names containing ISTD and add them to withheld columns
    istd <- grep("ISTD", metab_cols, value = TRUE)
    withheld_cols <- c(withheld_cols, istd)
    metab_cols <- setdiff(metab_cols, withheld_cols)
  }
  
  transformed_df <- df[, c(meta_cols, metab_cols)]
  
  if (transform == "none") {
    transform_str <- "None"
  } else if (transform == "log2") {
    transform_str <- "Log2"
    transformed_df[metab_cols] <- log(transformed_df[metab_cols], base = 2)
  } else if (transform == "TRN") {
    transform_str <- "TRN"
    transformed_df <- total_ratio_norm(transformed_df, metab_cols)
  }
  return(list(
    df = transformed_df,
    str = transform_str,
    withheld_cols = withheld_cols
  ))
}

group_stats <- function(df) {
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  names(df)[names(df) == "class"] <- "Group"
  df$batch <- NULL
  df$order <- NULL
  
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
  metab_cols <- setdiff(names(df), c("sample", "batch", "class", "order"))
  
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

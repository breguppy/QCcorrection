# plotting helpers
source("R/processing_helpers.R")
source("R/met_scatter_loess.R")
source("R/met_scatter_rf.R")
source("R/plotting_rsd_comparisons.R")
source("R/plotting_pca_comparisons.R")

# Helper for creating metabolite scatter plot
make_met_scatter <- function(d, met_col) {
  # choose the correct plotting function based on the correction method.
  cor_method <- d$corrected$str
  tryCatch({
    if (cor_method %in% c("Random Forest","Batchwise Random Forest")) {
      met_scatter_rf(d$filtered$df, d$filtered_corrected$df, i = met_col)
    } else if (cor_method %in% c("LOESS","Batchwise LOESS")) {
      met_scatter_loess(d$filtered$df, d$filtered_corrected$df, i = met_col)
    } else {
      ggplot2::ggplot() + ggplot2::labs(title = "No correction method selected.")
    }
  }, error = function(e) {
    showNotification(paste("Scatter failed:", e$message),
                     type = "error", duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "Scatter failed — see notification")
  })
}

# Helper for creating the RSD plot
make_rsd_plot <- function(p, d) {
  
  df_before <- d$filtered$df
  # Determine df_after based on rsd_compare selected by user.
  if (p$rsd_compare == "filtered_cor_data"){
    df_after <- d$filtered_corrected$df
    compared_to <- "Correction"
  } else {
    df_after <- d$transformed$df
    compared_to <- "Correction and Transformation"
  }
  
  # Need at least 1 metabolite column
  validate(
    need(ncol(df_before) > 4L, "No metabolites left before correction."),
    need(ncol(df_after)  > 4L, "No metabolites left after correction.")
  )
  
  tryCatch({
    if (identical(p$rsd_cal, "met")) {
      plot_rsd_comparison(df_before, df_after, compared_to)
    } else {
      plot_rsd_comparison_class_met(df_before, df_after, compared_to)
    }
    }, error = function(e) {
      showNotification(paste("RSD comparison failed:", e$message),
                       type = "error", duration = 8)
      ggplot2::ggplot() + ggplot2::labs(title = "RSD comparison failed — see notification")
    })
}

# Helper for creating the PCA plot
make_pca_plot <- function(p, d) {
  # get after based on pca_compare selected by user.
  if (p$pca_compare == "filtered_cor_data"){
    df <- d$filtered_corrected$df
    after <- d$filtered_corrected
    compared_to <- "Correction"
  } else {
    df <- d$transformed$df
    after <- d$transformed
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(df), c("sample","batch","class","order"))
  validate(
    need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    need(nrow(df) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- df[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  validate(need(any(keep), "All metabolite columns are constant/invalid after filtering."))
  
  # before data cannot have any missing values.
  before <- d$imputed
  tryCatch({
    plot_pca(p, before, after, compared_to)
  }, error = function(e) {
    showNotification(paste("PCA failed:", e$message),
                     type = "error", duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "PCA failed — see notification")
  })
}

# Helper function to ensure proper figure names
sanitize_figname <- function(name) {
  name <- gsub("[<>:\"/\\\\|?*]", "_", name)
  name <- gsub("\\s+", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_+|_+$", "", name)
  return(name)
}

# TODO: Cumulative RSD curves for other visualization. Still doesn't work
cumulative_met_rsd_auc <- function(df_after) {
  rsd_df <- metabolite_rsd(df_after)
  
  # Filter Missing value RSD
  rsd_sample <- rsd_df %>% filter(!is.na(RSD_NonQC))
  rsd_qc <- rsd_df %>% filter(!is.na(RSD_QC))
  
  # Compute ECDF for samples
  ecdf_sample <- rsd_sample %>%
    arrange(RSD_NonQC) %>%
    mutate(CDF = ecdf(RSD_NonQC)(RSD_NonQC))
  
  # Compute ECDF for samples
  ecdf_QC <- rsd_QC %>%
    arrange(RSD_QC) %>%
    mutate(CDF = ecdf(RSD_QC)(RSD_QC))
  # Compute mean RSD per batch and label
  summary_df <- rsd_df %>%
    group_by(Batch) %>%
    summarise(mean_RSD = mean(RSD, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_RSD) %>%
    mutate(
      Batch = as.character(Batch),
      BatchLabel = paste0(Batch, " (Mean RSD: ", round(mean_RSD, 1), "%)"),
      Batch = factor(Batch, levels = Batch)
    )
  
  # Merge labels into ECDF data
  ecdf_df <- ecdf_df %>%
    mutate(Batch = as.character(Batch)) %>%
    left_join(summary_df %>% select(Batch, BatchLabel), by = "Batch") %>%
    mutate(BatchLabel = factor(BatchLabel, levels = unique(summary_df$BatchLabel)))
  
  # Plot cumulative frequency curves
  p <- ggplot(ecdf_df, aes(x = RSD, y = CDF, color = BatchLabel)) +
    geom_line(size = 1) +
    labs(
      title = "Cumulative Frequency Curve of RSDs by Batch",
      x = "Relative Standard Deviation (RSD%)",
      y = "% of Metabolites",
      color = NULL
    ) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  today_str <- format(Sys.Date(), "%Y-%m-%d")
  fig_name <- paste0(today_str, "/Ref Batch/", cell_line, " ", data_type, " ", dose)
  ggsave(
    paste(fig_name, ".png"),
    plot = p,
    width = 8,
    height = 5,
    bg = "white"
  )
  ggsave(
    paste(fig_name, ".pdf"),
    plot = p,
    width = 8,
    height = 5
  )
  return(list(rsd = rsd_df, summary = summary_df))
}
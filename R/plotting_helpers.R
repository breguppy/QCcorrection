# plotting helpers
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tidyverse)
library(grid)
library(gridExtra)
source("R/processing_helpers.R")

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
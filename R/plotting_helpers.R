# ploting helpers

library(ggplot2)
library(patchwork)
library(dplyr)
source("processing_helpers.R")

met_scatter_rf <- function(data_raw, data_cor, i, metadata_cols = c("sample", "batch", "class", "order")) {
  
  qc_idx <- which(is.na(data_raw$class))
  nonqc_idx <- setdiff(seq_len(nrow(data_raw)), qc_idx)
  
  mean_before <- mean(data_raw[qc_idx, i])
  sd_before <- sd(data_raw[qc_idx, i])
  
  mean_after <- mean(data_cor[qc_idx, i])
  sd_after <- sd(data_cor[qc_idx, i])
  
  raw_batch_ranges <- data_raw %>%
    group_by(batch) %>%
    summarize(xmin = min(order), xmax = max(order), .groups = "drop") %>%
    arrange(xmin) %>%
    mutate(fill = rep(c("lightgray", "white"), length.out = n()))
  
  p_before <- ggplot() +
    # Background batch rectangles
    geom_rect(data = raw_batch_ranges,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.3) +
    scale_fill_identity() +
    # +-1 SD 
    geom_hline(yintercept = mean_before + sd_before,
               linetype = "dashed", colour = "grey20") +
    geom_hline(yintercept = mean_before - sd_before,
               linetype = "dashed", colour = "grey20") +
    # +-2 SD 
    geom_hline(yintercept = mean_before + 2*sd_before,
               linetype = "solid", colour = "#D55E00") +
    geom_hline(yintercept = mean_before - 2*sd_before,
               linetype = "solid", colour = "#D55E00") +
    # non-QC points in yellow
    geom_point(data = data_raw[nonqc_idx, ],
               aes(x = order, y = .data[[i]]),
               color = "#F5C710", size = 1.5) +
    # QC points in blue
    geom_point(data = data_raw[qc_idx, ],
               aes(x = order, y = .data[[i]]),
               color = "#305CDE", size = 2) +
    labs(
      title = "Raw",
      x = "Injection Order",
      y = "Intensity",
    ) +
    legend("top", c("Sample", "QC"), col = c("#F5C710", "#305CDE"), lty = 1, pch = 19, bty = "n", cex = 0.75, 
           horiz = TRUE)
    theme_minimal()
    
    cor_batch_ranges <- data_cor %>%
      group_by(batch) %>%
      summarize(xmin = min(order), xmax = max(order), .groups = "drop") %>%
      arrange(xmin) %>%
      mutate(fill = rep(c("lightgray", "white"), length.out = n()))
    
    p_after <- ggplot() +
      # Background batch rectangles
      geom_rect(data = cor_batch_ranges,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
                inherit.aes = FALSE, alpha = 0.3) +
      scale_fill_identity() +
      # +-1 SD 
      geom_hline(yintercept = mean_after + sd_after,
                 linetype = "dashed", colour = "grey20") +
      geom_hline(yintercept = mean_after - sd_after,
                 linetype = "dashed", colour = "grey20") +
      # +-2 SD 
      geom_hline(yintercept = mean_after + 2*sd_after,
                 linetype = "solid", colour = "#D55E00") +
      geom_hline(yintercept = mean_after - 2*sd_after,
                 linetype = "solid", colour = "#D55E00") +
      # non-QC points in yellow
      geom_point(data = data_cor[nonqc_idx, ],
                 aes(x = order, y = .data[[i]]),
                 color = "#F5C710", size = 1.5) +
      # QC points in blue
      geom_point(data = data_cor[qc_idx, ],
                 aes(x = order, y = .data[[i]]),
                 color = "#305CDE", size = 2) +
      labs(
        title = "Corrected",
        x = "Injection Order",
        y = "Intensity",
      ) +
      legend("top", c("Sample", "QC"), col = c("#F5C710", "#305CDE"), lty = 1, pch = 19, bty = "n", cex = 0.75, 
             horiz = TRUE)
    theme_minimal()
    
    p_before + p_after + plot_layout(ncol = 1, widths = 3, heights = 1)
}



plot_rsd_comparison <- function(df_before, df_after, corMethod) {
  rsdBefore <- metabolite_rsd(df_before)
  rsdAfter <- metabolite_rsd(df_after)
  
  # Merge data on Treatment and Metabolite
  df <- rsdBefore %>%
    rename(rsd_before = cv) %>%
    inner_join(rsdAfter %>% rename(rsd_after = cv),
               by = c("class", "Metabolite"))
  
  # Categorize changes in CV
  df <- df %>%
    mutate(change = case_when(
      rsd_after > rsd_before ~ "Increased",
      rsd_after < rsd_before ~ "Decreased",
      TRUE ~ "No Change"
    ))
  
  # Force all levels to be present
  df$change <- factor(df$change, levels = c("Increased", "No Change", "Decreased"))
  
  # Calculate percentages
  total <- nrow(df)
  perc <- df %>%
    count(change, .drop = FALSE) %>%
    complete(change = factor(c("Increased", "No Change", "Decreased"),
                             levels = c("Increased", "No Change", "Decreased")),
             fill = list(n = 0)) %>%
    mutate(percent = round(n / total * 100, 1))
  
  # Labels
  label_map <- setNames(
    paste0(c(
      paste0("Increased CV after ", corMethod, ": "),
      paste0("No change after ", corMethod, ": "),
      paste0("Decreased CV after ", corMethod, ": ")
    ), perc$percent, "%"),
    levels(df$change)
  )
  
  # Color mapping
  color_values <- c("Increased" = "darkred", "No Change" = "gray", "Decreased" = "darkgreen")
  
  # Dummy data with numeric NA to force legend appearance
  dummy_data <- data.frame(
    cv_before = as.numeric(NA),
    cv_after = as.numeric(NA),
    change = factor(c("Increased", "No Change", "Decreased"),
                    levels = c("Increased", "No Change", "Decreased"))
  )
  
  # Plot
  title_name <- paste0("Comparison of RSD Before and After ", corMethod)
  ggplot(df, aes(x = rsd_before, y = rsd_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point() +
    geom_point(data = dummy_data, aes(x = rsd_before, y = rsd_after, color = change), show.legend = TRUE) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map),
      labels = label_map,
      name = "RSD Change"
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = title_name
    ) +
    theme_minimal()
}
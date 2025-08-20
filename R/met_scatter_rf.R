# Matabolite scatter plot for random forest corrected data

library(ggplot2)
library(patchwork)
library(dplyr)

met_scatter_rf <- function(data_raw, data_cor, i) {
  
  raw_qc_idx <- which(data_raw$class == "QC")
  raw_nonqc_idx <- setdiff(seq_len(nrow(data_raw)), raw_qc_idx)
  
  cor_qc_idx <- which(data_cor$class == "QC")
  cor_nonqc_idx <- setdiff(seq_len(nrow(data_cor)), cor_qc_idx)
  
  raw_vals <- data_raw[[i]]
  qc_vals_raw <- raw_vals[ raw_qc_idx ]
  
  mean_before <- mean(qc_vals_raw, na.rm = TRUE)
  sd_before   <- sd(qc_vals_raw,   na.rm = TRUE)
  
  # same for data_cor
  cor_vals <- data_cor[[i]]
  qc_vals_cor  <- cor_vals[ cor_qc_idx ]
  
  mean_after <- mean(qc_vals_cor, na.rm = TRUE)
  sd_after   <- sd(qc_vals_cor,   na.rm = TRUE)
  
  sd_df_before <- tibble::tibble(
    y  = c(mean_before +  sd_before,
           mean_before -  sd_before,
           mean_before + 2*sd_before,
           mean_before - 2*sd_before),
    sd = factor(c("±1 SD","±1 SD","±2 SD","±2 SD"),
                levels = c("±1 SD","±2 SD"))
  )
  
  sd_df_after <- tibble::tibble(
    y  = c(mean_after  +  sd_after,
           mean_after  -  sd_after,
           mean_after  + 2*sd_after,
           mean_after  - 2*sd_after),
    sd = factor(c("±1 SD","±1 SD","±2 SD","±2 SD"),
                levels = c("±1 SD","±2 SD"))
  )
  
  raw_batch_ranges <- data_raw %>%
    group_by(batch) %>%
    summarize(xmin = min(order), xmax = max(order), .groups = "drop") %>%
    arrange(xmin) %>%
    mutate(fill = rep(c("lightgray", "white"), length.out = n()))
  
  cor_batch_ranges <- data_cor %>%
    group_by(batch) %>%
    summarize(xmin = min(order), xmax = max(order), .groups = "drop") %>%
    arrange(xmin) %>%
    mutate(fill = rep(c("lightgray", "white"), length.out = n()))
  
  data_raw$type <- ifelse(data_raw$class == "QC", "QC", "Sample")
  data_cor$type <- ifelse(data_cor$class == "QC", "QC", "Sample")
  
  color_scale <- scale_color_manual(
    name   = "Type: ",
    values = c(Sample = "#F5C710", QC = "#305CDE")
  )
  lty_scale <- scale_linetype_manual(
    name   = "SD Range: ",
    values = c("±1 SD" = "dashed", "±2 SD" = "solid"),
    guide  = guide_legend(override.aes = list(color = c("grey20","#950606")))
  )
  
  p_before <- ggplot(data_raw, aes(x = order, y = .data[[i]], color = type), alpha = 0.8) +
    geom_rect(data = raw_batch_ranges,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.3) +
    scale_fill_identity() +
    geom_hline(data = sd_df_before, aes(yintercept = y, linetype = sd),
               color = ifelse(sd_df_before$sd=="±1 SD","grey20","#950606"), linewidth = 1) +
    geom_point(data = data_raw %>% filter(type == "Sample"),
               aes(order, .data[[i]], color = type), size = 2) +
    geom_point(data = data_raw %>% filter(type == "QC"),
               aes(order, .data[[i]], color = type), size = 2) +
    color_scale + lty_scale +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 10),
      legend.text= element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    labs(title = "Raw", x = "Injection Order", y = "Intensity")
  
  p_after <- ggplot(data_cor, aes(x = order, y = .data[[i]], color = type)) +
    geom_rect(data = cor_batch_ranges,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.3) +
    scale_fill_identity() +
    geom_hline(data = sd_df_after, aes(yintercept = y, linetype = sd),
               color = ifelse(sd_df_after$sd=="±1 SD","grey20","#950606"), linewidth = 1) +
    geom_point(data = data_cor %>% filter(type == "Sample"),
               aes(order, .data[[i]], color = type), size = 2) +
    geom_point(data = data_cor %>% filter(type == "QC"),
               aes(order, .data[[i]], color = type), size = 2) +
    color_scale + lty_scale +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 10),
      legend.text= element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    labs(title = "Corrected", x = "Injection Order", y = "Intensity")
  
  
  (p_before + p_after) +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(title = i) &
    theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "bottom")
}
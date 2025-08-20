# metabolite scatter plots for polynomial fit

library(ggplot2)
library(patchwork)
library(dplyr)

met_scatter_loess <- function(data_raw, data_cor, i) {
  
  raw_qc_idx <- which(data_raw$class == "QC")
  raw_nonqc_idx <- setdiff(seq_len(nrow(data_raw)), raw_qc_idx)
  
  cor_qc_idx <- which(data_cor$class == "QC")
  cor_nonqc_idx <- setdiff(seq_len(nrow(data_cor)), cor_qc_idx)
  
  raw_vals <- data_raw[[i]]
  qc_vals_raw <- raw_vals[ raw_qc_idx ]
  
  # same for data_cor
  cor_vals <- data_cor[[i]]
  qc_vals_cor  <- cor_vals[ cor_qc_idx ]
  
  
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
  
  p_before <- ggplot(data_raw, aes(x = order, y = .data[[i]], color = type), alpha = 0.8) +
    geom_rect(data = raw_batch_ranges,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              inherit.aes = FALSE, alpha = 0.3) +
    scale_fill_identity() +
    geom_point(data = data_raw %>% filter(type == "Sample"),
               aes(order, .data[[i]], color = type), size = 2) +
    geom_smooth(data = data_raw %>% filter(type=="QC"),
                aes(x = order, y = .data[[i]]), formula = "y ~ x",
                method = "loess", span  = 0.75,
                fill = "#305CDE", show.legend = FALSE) +
    geom_point(data = data_raw %>% filter(type == "QC"),
               aes(order, .data[[i]], color = type), size = 2) +
    color_scale +
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
    geom_point(data = data_cor %>% filter(type == "Sample"),
               aes(order, .data[[i]], color = type), size = 2) +
    geom_point(data = data_cor %>% filter(type == "QC"),
               aes(order, .data[[i]], color = type), size = 2) +
    color_scale +
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
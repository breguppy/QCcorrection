# ploting helpers

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
source("R/processing_helpers.R")


plot_rsd_comparison <- function(df_before, df_after) {
  rsdBefore <- metabolite_rsd(df_before)
  rsdAfter <- metabolite_rsd(df_after)
  
  # Merge data on Treatment and Metabolite
  df <- rsdBefore %>%
    rename(rsd_qc_before = RSD_QC, rsd_nonqc_before = RSD_NonQC) %>%
    inner_join(rsdAfter %>% rename(rsd_qc_after = RSD_QC, rsd_nonqc_after = RSD_NonQC),
               by = c("Metabolite"))
  
  # Categorize changes in CV
  df_samples <- df %>%
    mutate(change = case_when(
      rsd_nonqc_after > rsd_nonqc_before ~ "Increased",
      rsd_nonqc_after < rsd_nonqc_before ~ "Decreased",
      TRUE ~ "No Change"
    ))
  
  df_qcs <- df %>%
    mutate(change = case_when(
      rsd_qc_after > rsd_qc_before ~ "Increased",
      rsd_qc_after < rsd_qc_before ~ "Decreased",
      TRUE ~ "No Change"
    ))
  # Force all levels to be present
  df_samples$change <- factor(df_samples$change, levels = c("Increased", "No Change", "Decreased"))
  df_qcs$change <- factor(df_qcs$change, levels = c("Increased", "No Change", "Decreased"))
  
  # Calculate percentages
  total_samples <- nrow(df_samples)
  perc_samples <- df_samples %>%
    count(change, .drop = FALSE) %>%
    complete(change = factor(c("Increased", "No Change", "Decreased"),
                             levels = c("Increased", "No Change", "Decreased")),
             fill = list(n = 0)) %>%
    mutate(percent = round(n / total_samples * 100, 1))
  # Calculate percentages
  total_qcs <- nrow(df_qcs)
  perc_qcs <- df_qcs %>%
    count(change, .drop = FALSE) %>%
    complete(change = factor(c("Increased", "No Change", "Decreased"),
                             levels = c("Increased", "No Change", "Decreased")),
             fill = list(n = 0)) %>%
    mutate(percent = round(n / total_qcs * 100, 1))
  
  # Labels
  label_map_samples <- setNames(
    paste0(c(
      "Increased RSD: ",
      "No change: ",
      "Decreased RSD: "
    ), perc_samples$percent, "%"),
    levels(df_samples$change)
  )
  label_map_qcs <- setNames(
    paste0(c(
      "Increased RSD: ",
      "No change: ",
      "Decreased RSD: "
    ), perc_qcs$percent, "%"),
    levels(df_qcs$change)
  )
  
  # Color mapping
  color_values <- c("Increased" = "#B22222", "No Change" = "gray25", "Decreased" = "#234F1E")
  
  # Dummy data with numeric NA to force legend appearance
  dummy_data <- data.frame(
    rsd_qc_before = as.numeric(NA),
    rsd_qc_after = as.numeric(NA),
    rsd_nonqc_before = as.numeric(NA),
    rsd_nonqc_after = as.numeric(NA),
    change = factor(c("Increased", "No Change", "Decreased"),
                    levels = c("Increased", "No Change", "Decreased"))
  )
  
  # Plot
  title_name <- "Comparison of RSD Before and After Correction"
  p1 <- ggplot(df_samples, aes(x = rsd_nonqc_before, y = rsd_nonqc_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3) +
    geom_point(data = dummy_data, 
               aes(x = rsd_nonqc_before, y = rsd_nonqc_after, color = change), 
               show.legend = TRUE, size = 3) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map_samples),
      labels = label_map_samples,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = c(0.24, 0.89),
      legend.background = element_rect(fill="white", 
                                       size=0.5, linetype="solid"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "Samples"
    )
  
  p2 <- ggplot(df_qcs, aes(x = rsd_qc_before, y = rsd_qc_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3) +
    geom_point(data = dummy_data, 
               aes(x = rsd_qc_before, y = rsd_qc_after, color = change), 
               show.legend = TRUE, size = 3) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map_qcs),
      labels = label_map_qcs,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = c(0.24, 0.89),
      legend.background = element_rect(fill="white", 
                                       size=0.5, linetype="solid"),
      panel.border = element_rect(color = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "QCs"
    )
  
  p1 + p2 + plot_annotation(title_name, theme=theme(plot.title=element_text(size = 30, face = "bold", hjust=0.5)))
}

plot_rsd_comparison_class_met <- function(df_before, df_after) {
  rsdBefore <- class_metabolite_rsd(df_before)
  rsdAfter <- class_metabolite_rsd(df_after)
  
  # Merge data on Treatment and Metabolite
  df <- rsdBefore %>%
    rename(rsd_before = RSD) %>%
    inner_join(rsdAfter %>% rename(rsd_after = RSD),
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
  df_qcs <- df[df$class == "QC", ]
  df_samples <- df[df$class != "QC", ]
  
  total_samples <- nrow(df_samples)
  perc_samples <- df_samples %>%
    count(change, .drop = FALSE) %>%
    complete(change = factor(c("Increased", "No Change", "Decreased"),
                             levels = c("Increased", "No Change", "Decreased")),
             fill = list(n = 0)) %>%
    mutate(percent = round(n / total_samples * 100, 1))
  # Calculate percentages
  total_qcs <- nrow(df_qcs)
  perc_qcs <- df_qcs %>%
    count(change, .drop = FALSE) %>%
    complete(change = factor(c("Increased", "No Change", "Decreased"),
                             levels = c("Increased", "No Change", "Decreased")),
             fill = list(n = 0)) %>%
    mutate(percent = round(n / total_qcs * 100, 1))
  
  # Labels
  label_map_samples <- setNames(
    paste0(c(
      "Increased RSD: ",
      "No change: ",
      "Decreased RSD: "
    ), perc_samples$percent, "%"),
    levels(df_samples$change)
  )
  label_map_qcs <- setNames(
    paste0(c(
      "Increased RSD: ",
      "No change: ",
      "Decreased RSD: "
    ), perc_qcs$percent, "%"),
    levels(df_qcs$change)
  )
  
  # Color mapping
  color_values <- c("Increased" = "#B22222", "No Change" = "gray25", "Decreased" = "#234F1E")
  
  # Dummy data with numeric NA to force legend appearance
  dummy_data <- data.frame(
    rsd_before = as.numeric(NA),
    rsd_after = as.numeric(NA),
    change = factor(c("Increased", "No Change", "Decreased"),
                    levels = c("Increased", "No Change", "Decreased"))
  )
  
  # Plot
  title_name <- "Comparison of RSD Before and After Correction"
  p1 <- ggplot(df_samples, aes(x = rsd_before, y = rsd_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3) +
    geom_point(data = dummy_data, 
               aes(x = rsd_before, y = rsd_after, color = change), 
               show.legend = TRUE, size = 3) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map_samples),
      labels = label_map_samples,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = c(0.24, 0.89),
      legend.background = element_rect(fill="white", 
                                       size=0.5, linetype="solid"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "Samples"
    )
  
  p2 <- ggplot(df_qcs, aes(x = rsd_before, y = rsd_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3) +
    geom_point(data = dummy_data, 
               aes(x = rsd_before, y = rsd_after, color = change), 
               show.legend = TRUE, size = 3) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map_qcs),
      labels = label_map_qcs,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = c(0.24, 0.89),
      legend.background = element_rect(fill="white", 
                                       size=0.5, linetype="solid"),
      panel.border = element_rect(color = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "QCs"
    )
  
  p1 + p2 + plot_annotation(title_name, theme=theme(plot.title=element_text(size = 30, face = "bold", hjust=0.5)))
}
# ploting helpers

library(ggplot2)
library(patchwork)
library(dplyr)
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
  df <- df %>%
    mutate(change = case_when(
      rsd_qc_after > rsd_qc_before ~ "Increased",
      rsd_qc_after < rsd_qc_before ~ "Decreased",
      rsd_nonqc_after > rsd_nonqc_before ~ "Increased",
      rsd_nonqc_after < rsd_nonqc_before ~ "Decreased",
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
      "Increased CV after correction: ",
      "No change after correction: ",
      "Decreased CV after correction: "
    ), perc$percent, "%"),
    levels(df$change)
  )
  
  # Color mapping
  color_values <- c("Increased" = "darkred", "No Change" = "gray", "Decreased" = "darkgreen")
  
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
  p1 <- ggplot(df, aes(x = rsd_nonqc_before, y = rsd_nonqc_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point() +
    geom_point(data = dummy_data, 
               aes(x = rsd_nonqc_before, y = rsd_nonqc_after, color = change), 
               show.legend = TRUE) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map),
      labels = label_map,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "Samples"
    )
  
  p2 <- ggplot(df, aes(x = rsd_qc_before, y = rsd_qc_after, color = change)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point() +
    geom_point(data = dummy_data, 
               aes(x = rsd_qc_before, y = rsd_qc_after, color = change), 
               show.legend = TRUE) +
    scale_color_manual(
      values = color_values,
      breaks = names(label_map),
      labels = label_map,
      name = "RSD Change"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.text= element_text(size = 16),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill=NA, linewidth=1)
    ) +
    labs(
      x = "RSD Before",
      y = "RSD After",
      title = "QCs"
    )
  
  p1 + p2 + plot_annotation(title_name, theme=theme(plot.title=element_text(hjust=0.5)))
}
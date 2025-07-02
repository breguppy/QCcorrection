# ploting helpers

library(ggplot2)
library(patchwork)
library(dplyr)
R_files <- list.files(path = getwd(), pattern = "\\.R$", full.names = TRUE)
sapply(R_files, source)


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
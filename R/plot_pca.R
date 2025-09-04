# plotting PCA:

# PCA plots coloring by color_col
plot_pca <- function(p, before, after, compared_to) {
  # Get overlapping metabolite columns
  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols <- intersect(setdiff(names(before$df), meta_cols), setdiff(names(after$df), meta_cols))
  
  # impute missing values in PCA if there are any
  if (any(is.na(after$df[metab_cols]))) {
    results <- impute_missing(after$df, metab_cols, p$qcImputeM, p$samImputeM)
    after <- results$df
  } else {
    after <- after$df
  }
  
  # before PCA
  before_pca_result <- prcomp(before$df[, metab_cols], center = TRUE, scale. = TRUE)
  before_pca_df <- as.data.frame(before_pca_result$x[, 1:2])
  colnames(before_pca_df) <- c("PC1", "PC2")
  before_pca_df <- bind_cols(before_pca_df, before$df[meta_cols]) #%>% mutate(panel = "Before")
  
  # After PCA
  after_pca_result <- prcomp(after[, metab_cols], center = TRUE, scale. = TRUE)
  after_pca_df <- as.data.frame(after_pca_result$x[, 1:2])
  colnames(after_pca_df) <- c("PC1", "PC2")
  after_pca_df <- bind_cols(after_pca_df, after[meta_cols]) #%>% mutate(panel = "After")
  
  # Combine scores for consistent axis scaling
  combined <- bind_rows(before_pca_df, after_pca_df)
  x_limits <- range(combined$PC1)
  y_limits <- range(combined$PC2)
  
  # Variance explained
  var_exp_raw <- summary(before_pca_result)$importance[2, 1:2] * 100
  var_exp_cor <- summary(after_pca_result)$importance[2, 1:2] * 100
  ann_df <- data.frame(
    panel = c("Before", "Before", "After", "After"),
    label = c(
      sprintf("PC1 (%.1f%%)", var_exp_raw[1]),
      sprintf("PC2 (%.1f%%)", var_exp_raw[2]),
      sprintf("PC1 (%.1f%%)", var_exp_cor[1]),
      sprintf("PC2 (%.1f%%)", var_exp_cor[2])
    ),
    x = c(
      mean(x_limits),
      x_limits[1] - diff(x_limits) * 0.15,
      mean(x_limits),
      x_limits[1] - diff(x_limits) * 0.15
    ),
    y = c(
      y_limits[1] - diff(y_limits) * 0.10,
      mean(y_limits),
      y_limits[1] - diff(y_limits) * 0.10,
      mean(y_limits)
    ),
    hjust = c(0.5, 0),
    vjust = c(1, 0.5)
  )
  
  # Custom color palette
  cbPalette <- c(
    "#F3C300",
    "#875692",
    "#ee7733",
    "#A1CAF1",
    "#BE0032",
    "#C2B280",
    "#555555",
    "#008856",
    "#E68FAC",
    "#0067A5",
    "#F99379",
    "#332288",
    "#F6A600",
    "#B3446C",
    "#DCD300",
    "#882D17",
    "#8DB600",
    "#654522",
    "#E25822",
    "#2B3D26",
    "#bbbbbb",
    "#000000",
    "#33bbee",
    "#ccddaa",
    "#225555"
  )
  
  color_levels <- sort(unique(c(before_pca_df[[p$color_col]], after_pca_df[[p$color_col]])))
  if (length(color_levels) > length(cbPalette)) {
    stop("Not enough colors in cbPalette for the number of groups in color_col.")
  }
  color_values <- setNames(cbPalette[seq_along(color_levels)], color_levels)
  combined[[p$color_col]] <- factor(combined[[p$color_col]], levels = color_levels)
  before_pca_df[[p$color_col]] <- factor(before_pca_df[[p$color_col]], levels = color_levels)
  after_pca_df [[p$color_col]] <- factor(after_pca_df [[p$color_col]], levels = color_levels)
  
  
  # Theme with larger fonts
  big_font_theme <- theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  p1 <- ggplot(before_pca_df, aes(x = PC1, y = PC2, color = .data[[p$color_col]])) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = "Before",
      x = paste0("PC1 (", round(var_exp_raw[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp_raw[2], 1), "%)")
    ) +
    xlim(x_limits) + ylim(y_limits) +
    scale_color_manual(
      values = color_values,
      name = p$color_col,
      drop = FALSE,
      na.translate = FALSE
    ) +
    big_font_theme +
    theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  # Corrected plot with no legend
  p2 <- ggplot(after_pca_df, aes(x = PC1, y = PC2, color = .data[[p$color_col]])) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = "After",
      x = paste0("PC1 (", round(var_exp_cor[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp_cor[2], 1), "%)")
    ) +
    xlim(x_limits) + ylim(y_limits) +
    scale_color_manual(values = color_values,
                       drop = FALSE,
                       na.translate = FALSE) +
    big_font_theme +
    theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  p_leg <- ggplot(before_pca_df, aes(PC1, PC2, color = .data[[p$color_col]])) +
    geom_point(size = 2) +
    scale_color_manual(
      values = color_values,
      name = p$color_col,
      drop = FALSE,
      na.translate = FALSE
    ) +
    guides(
      fill = "none",
      size = "none",
      shape = "none",
      alpha = "none",
      linetype = "none"
    ) +
    big_font_theme +
    theme(legend.position = "right",
          legend.box.margin = ggplot2::margin(0, 0, 0, 0))
  
  legend <- cowplot::get_legend(p_leg)
  
  p_combined <- cowplot::plot_grid(
    p1,
    p2,
    legend,
    labels = NULL,
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(1, 1, 0.22)
  )
  
  p_with_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste("Comparison of PCA Before and After", compared_to),
      fontface = "bold",
      x = 0.5,
      y = 0.98,
      hjust = 0.5,
      vjust = 1,
      size = 14
    ) +
    cowplot::draw_plot(
      p_combined,
      x = 0,
      y = 0,
      width = 1,
      height = 0.93
    )
  
  return(p_with_title)
}
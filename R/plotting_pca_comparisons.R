# plotting PCA:

# PCA plots coloring by color_col
plot_pca <- function(input,
                     before,
                     after,
                     color_col) {
  meta_cols <- c("sample", "batch", "class", "order")
  metab_cols1 <- setdiff(names(before$df), meta_cols)
  metab_cols2 <- setdiff(names(after$df), meta_cols)
  metab_cols <- intersect(metab_cols1, metab_cols2)
  
  if (any(is.na(after$df[metab_cols]))) {
    results <- impute_missing(after$df,
                              metab_cols,
                              input$qcImputeM,
                              input$samImputeM)
    after <- results$df
  } else {
    after <- after$df
  }
  # before PCA
  before_pca_result <- prcomp(before$df[, metab_cols], center = TRUE, scale. = TRUE)
  before_pca_df <- as.data.frame(before_pca_result$x[, 1:2])
  colnames(before_pca_df) <- c("PC1", "PC2")
  before_pca_df <- bind_cols(before_pca_df, before$df[meta_cols])
  
  # After PCA
  after_pca_result <- prcomp(after[, metab_cols], center = TRUE, scale. = TRUE)
  after_pca_df <- as.data.frame(after_pca_result$x[, 1:2])
  colnames(after_pca_df) <- c("PC1", "PC2")
  after_pca_df <- bind_cols(after_pca_df, after[meta_cols])
  
  # Combine scores for consistent axis scaling
  combined <- bind_rows(before_pca_df, after_pca_df)
  x_limits <- range(combined$PC1)
  y_limits <- range(combined$PC2)
  
  # Variance explained
  var_exp_raw <- summary(before_pca_result)$importance[2, 1:2] * 100
  var_exp_cor <- summary(after_pca_result)$importance[2, 1:2] * 100
  
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
  
  color_levels <- unique(c(before_pca_df[[color_col]], after_pca_df[[color_col]]))
  color_levels <- sort(color_levels)
  n_levels <- length(color_levels)
  if (n_levels > length(cbPalette)) {
    stop("Not enough colors in cbPalette for the number of groups in color_col.")
  }
  color_values <- setNames(cbPalette[1:n_levels], color_levels)
  before_pca_df[[color_col]] <- factor(before_pca_df[[color_col]], levels = color_levels)
  after_pca_df[[color_col]] <- factor(after_pca_df[[color_col]], levels = color_levels)
  
  # Theme with larger fonts
  big_font_theme <- theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13)
    )
  
  # Raw plot
  p1 <- ggplot(before_pca_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 4, alpha = 0.8) +
    labs(
      title = paste0("PCA - Before"),
      x = paste0("PC1 (", round(var_exp_raw[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp_raw[2], 1), "%)")
    ) +
    xlim(x_limits) + ylim(y_limits) +
    scale_color_manual(values = cbPalette, name = color_col) +
    big_font_theme +
    theme(legend.position = "right")
  
  # Corrected plot with no legend
  p2 <- ggplot(after_pca_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 4, alpha = 0.8) +
    labs(
      title = paste0("PCA - After"),
      x = paste0("PC1 (", round(var_exp_cor[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp_cor[2], 1), "%)")
    ) +
    xlim(x_limits) + ylim(y_limits) +
    scale_color_manual(values = cbPalette) +
    big_font_theme +
    theme(legend.position = "none")
  
  p <- p1 + p2 + plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  return(p)
}
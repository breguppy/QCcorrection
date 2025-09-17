#' PCA plotting
#' @keywords internal
#' @noRd
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
  before_p <- stats::prcomp(before$df[, metab_cols], center = TRUE, scale. = TRUE)
  before_df <- tibble::as_tibble(before_p$x[, 1:2, drop = FALSE], .name_repair = "minimal")
  names(before_df) <- c("PC1", "PC2")
  before_df <- dplyr::bind_cols(before_df, before$df[meta_cols])
  
  # After PCA
  after_p  <- stats::prcomp(after[, metab_cols], center = TRUE, scale. = TRUE)
  after_df  <- tibble::as_tibble(after_p$x[, 1:2, drop = FALSE], .name_repair = "minimal")
  names(after_df) <- c("PC1", "PC2")
  after_df  <- dplyr::bind_cols(after_df, after[meta_cols])
  
  # Combine scores for consistent axis scaling
  combined <- dplyr::bind_rows(before_df, after_df)
  x_limits <- range(combined$PC1)
  y_limits <- range(combined$PC2)
  
  # Variance explained
  var_raw <- summary(before_p)$importance[2, 1:2] * 100
  var_cor <- summary(after_p)$importance[2, 1:2] * 100
  
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
  
  col <- p$color_col %||% "class"
  lvls <- sort(unique(c(before_df[[col]], after_df[[col]])))
  if (length(lvls) > length(cbPalette))
    stop("Too many groups for palette.")
  cols <- stats::setNames(cbPalette[seq_along(lvls)], lvls)
  combined[[col]]  <- factor(combined[[col]], levels = lvls)
  before_df[[col]] <- factor(before_df[[col]], levels = lvls)
  after_df[[col]]  <- factor(after_df[[col]], levels = lvls)
  
  # Theme with larger fonts
  big_font_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 14,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text  = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10)
    )
  
  # before plot
  p1 <- ggplot2::ggplot(before_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "Before",
      x = sprintf("PC1 (%.1f%%)", var_raw[1]),
      y = sprintf("PC2 (%.1f%%)", var_raw[2])
    ) +
    ggplot2::xlim(x_limits) + ggplot2::ylim(y_limits) +
    ggplot2::scale_color_manual(
      values = cols,
      name = col,
      drop = FALSE,
      na.translate = FALSE
    ) +
    big_font_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  # after plot
  p2 <- ggplot2::ggplot(after_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "After",
      x = sprintf("PC1 (%.1f%%)", var_cor[1]),
      y = sprintf("PC2 (%.1f%%)", var_cor[2])
    ) +
    ggplot2::xlim(x_limits) + ggplot2::ylim(y_limits) +
    ggplot2::scale_color_manual(values = cols,
                                drop = FALSE,
                                na.translate = FALSE) +
    big_font_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      plot.margin = ggplot2::margin(10, 5, 10, 5)
    )
  
  # Legend plot
  p_leg <- ggplot2::ggplot(before_df, ggplot2::aes(PC1, PC2, color = .data[[col]])) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(
      values = cols,
      name = col,
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::guides(
      fill = "none",
      size = "none",
      shape = "none",
      alpha = "none",
      linetype = "none"
    ) +
    big_font_theme + ggplot2::theme(legend.position = "right",
                                    legend.box.margin = ggplot2::margin(0, 0, 0, 0))
  
  leg <- cowplot::get_plot_component(p_leg, "guide-box", return_all = TRUE)[[1]]
  comb <- cowplot::plot_grid(
    p1,
    p2,
    leg,
    nrow = 1,
    rel_widths = c(1, 1, 0.3),
    labels = NULL,
    align = "hv",
    axis = "tblr"
  )
  
  cowplot::ggdraw() +
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
      comb,
      x = 0,
      y = 0,
      width = 1,
      height = 0.93
    )
}
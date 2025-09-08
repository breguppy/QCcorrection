#' Metabolite scatter plot for random forest corrected data
#' @keywords internal
#' @noRd

met_scatter_rf <- function(data_raw, data_cor, i) {
  data_raw <- dplyr::mutate(
    data_raw,
    type  = ifelse(class == "QC", "QC", "Sample"),
    panel = factor("Raw", levels = c("Raw","Corrected"))
  )
  data_cor <- dplyr::mutate(
    data_cor,
    type  = ifelse(class == "QC", "QC", "Sample"),
    panel = factor("Corrected", levels = c("Raw","Corrected"))
  )
  df_all <- dplyr::bind_rows(data_raw, data_cor)
  
  # ensure numeric order
  df_all$order   <- suppressWarnings(as.numeric(df_all$order))
  data_raw$order <- suppressWarnings(as.numeric(data_raw$order))
  data_cor$order <- suppressWarnings(as.numeric(data_cor$order))
  
  # QC SD bands (guard for empty/NA)
  get_stats <- function(df, panel) {
    qc <- df[df$type == "QC", i, drop = TRUE]
    m  <- mean(qc, na.rm = TRUE)
    s  <- stats::sd(qc, na.rm = TRUE)
    if (!is.finite(m) || !is.finite(s)) return(tibble::tibble(y = numeric(0), sd = factor(), panel = factor()))
    tibble::tibble(
      y     = c(m + s, m - s, m + 2*s, m - 2*s),
      sd    = factor(c("±1 SD","±1 SD","±2 SD","±2 SD"), levels = c("±1 SD","±2 SD")),
      panel = factor(panel, levels = c("Raw","Corrected"))
    )
  }
  sd_df <- dplyr::bind_rows(get_stats(data_raw, "Raw"),
                            get_stats(data_cor, "Corrected"))
  
  # batch shading with safeguards
  get_batches <- function(df, panel) {
    rng <- df |>
      dplyr::filter(!is.na(order)) |>
      dplyr::group_by(batch) |>
      dplyr::summarise(
        xmin = suppressWarnings(min(as.numeric(order), na.rm = TRUE)),
        xmax = suppressWarnings(max(as.numeric(order), na.rm = TRUE)),
        .groups = "drop"
      ) |>
      dplyr::filter(is.finite(xmin), is.finite(xmax), xmax >= xmin)
    
    if (nrow(rng) < 2L) return(rng[0, ])  # skip shading for 0/1 batch
    
    rng |>
      dplyr::arrange(xmin) |>
      dplyr::mutate(
        fill  = rep(c("lightgray","white"), length.out = dplyr::n()),
        panel = factor(panel, levels = c("Raw","Corrected"))
      )
  }
  
  batch_ranges <- dplyr::bind_rows(
    get_batches(data_raw, "Raw"),
    get_batches(data_cor, "Corrected")
  )
  
  color_scale <- ggplot2::scale_color_manual(
    name = "Type:", values = c(Sample = "#F5C710", QC = "#305CDE")
  )
  lty_scale <- ggplot2::scale_linetype_manual(
    name = "SD Range:",
    values = c("±1 SD" = "dashed", "±2 SD" = "solid"),
    guide  = ggplot2::guide_legend(override.aes = list(color = c("grey20", "#950606")))
  )
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = order, y = .data[[i]]))
  
  if (nrow(batch_ranges)) {
    p <- p +
      ggplot2::geom_rect(
        data = batch_ranges,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax, fill = fill),
        ymin = -Inf, ymax = Inf,
        inherit.aes = FALSE, alpha = 0.3, show.legend = FALSE
      ) +
      ggplot2::scale_fill_identity(guide = "none")
  }
  
  # SD lines only if we have rows
  if (nrow(sd_df)) {
    p <- p +
      ggplot2::geom_hline(
        data = sd_df |> dplyr::filter(sd == "±1 SD"),
        ggplot2::aes(yintercept = y, linetype = sd),
        color = "grey20", linewidth = 0.75
      ) +
      ggplot2::geom_hline(
        data = sd_df |> dplyr::filter(sd == "±2 SD"),
        ggplot2::aes(yintercept = y, linetype = sd),
        color = "#950606", linewidth = 0.75
      )
  }
  
  p +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "Sample"),
      ggplot2::aes(color = type), size = 2, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "QC"),
      ggplot2::aes(color = type), size = 2, na.rm = TRUE
    ) +
    color_scale + lty_scale +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(title = i, x = "Injection Order", y = "Intensity") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 12),
      axis.text    = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      legend.box        = "horizontal",
      legend.box.just   = "center",
      legend.key.width  = grid::unit(0.5, "cm"),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.margin     = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      legend.box.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.5, 0.5, 0.8, 0.5), "cm")
    )
}

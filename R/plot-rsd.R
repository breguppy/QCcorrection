#' RSD comparison plots
#' @keywords internal
#' @noRd
lab_levels   <- c("Increased", "No Change", "Decreased")

#' @keywords internal
#' @noRd
color_values <- c(
  "Increased" = "#B22222",
  "No Change" = "gray25",
  "Decreased" = "#234F1E"
)

#' @keywords internal
#' @noRd
# RSD comparison scatter plots when computing RSD by metabolite
plot_rsd_comparison <- function(df_before, df_after, compared_to) {
  rsdBefore <- metabolite_rsd(df_before)
  rsdAfter <- metabolite_rsd(df_after)
  
  # Merge data on Metabolite
  df <- dplyr::inner_join(
    dplyr::rename(
      rsdBefore,
      rsd_qc_before = RSD_QC,
      rsd_nonqc_before = RSD_NonQC
    ),
    dplyr::rename(
      rsdAfter,
      rsd_qc_after  = RSD_QC,
      rsd_nonqc_after  = RSD_NonQC
    ),
    by = "Metabolite"
  )
  # If empty, return empty plot
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No overlapping metabolites"))
  }
  
  # Get only finite data and categorize changes
  df_samples <- df %>%
    dplyr::filter(is.finite(rsd_nonqc_before), is.finite(rsd_nonqc_after)) %>%
    tag_changes("rsd_nonqc_before", "rsd_nonqc_after") %>%
    dplyr::transmute(Type = "Samples",
                     before = rsd_nonqc_before,
                     after = rsd_nonqc_after,
                     change)
  
  df_qcs <- df %>%
    dplyr::filter(is.finite(rsd_qc_before), is.finite(rsd_qc_after)) %>%
    tag_changes("rsd_qc_before", "rsd_qc_after") %>%
    dplyr::transmute(Type = "QCs",
                     before = rsd_qc_before,
                     after = rsd_qc_after,
                     change)
  
  d_all <- dplyr::bind_rows(df_samples, df_qcs) %>%
    dplyr::mutate(Type = factor(Type, levels = c("Samples", "QCs")),
                  change = factor(change, levels = lab_levels))
  
  # make legend map
  facet_labs <- facet_label_map(d_all)
  
  # single faceted ggplot
  mk_plot(d_all, "before", "after", facet_labs, compared_to)
}

#' @keywords internal
#' @noRd
# RSD comparison scatter plot when computing RSD by metabolite-class
plot_rsd_comparison_class_met <- function(df_before, df_after, compared_to) {
  rsdBefore <- class_metabolite_rsd(df_before)
  rsdAfter <- class_metabolite_rsd(df_after)
  
  # Merge data on class and Metabolite
  df <- dplyr::inner_join(
    dplyr::rename(rsdBefore, rsd_before = RSD),
    dplyr::rename(rsdAfter, rsd_after  = RSD),
    by = c("class", "Metabolite")
  )
  # If empty, return empty plot
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No overlapping metabolites"))
  }
  
  # Get only finite data and Categorize changes
  df_samples <- df[df$class != "QC", ] %>% dplyr::filter(is.finite(rsd_before), is.finite(rsd_after)) %>%
    tag_changes("rsd_before", "rsd_after") %>% dplyr::transmute(Type = "Samples",
                                                                before = rsd_before,
                                                                after = rsd_after,
                                                                change)
  
  df_qcs <- df[df$class == "QC", ] %>% dplyr::filter(is.finite(rsd_before), is.finite(rsd_after)) %>%
    tag_changes("rsd_before", "rsd_after") %>% dplyr::transmute(Type = "QCs",
                                                                before = rsd_before,
                                                                after = rsd_after,
                                                                change)
  
  d_all <- dplyr::bind_rows(df_samples, df_qcs) %>%
    dplyr::mutate(Type = factor(Type, levels = c("Samples", "QCs")),
                  change = factor(change, levels = lab_levels))
  
  # make legend map
  facet_labs <- facet_label_map(d_all)
  
  # single faceted ggplot
  mk_plot(d_all, "before", "after", facet_labs, compared_to)
}


#' Plot RSD distributions before and after correction
#'
#' @param df_before Data frame for computing Raw RSD
#' @param df_after Data frame for computing met RSD after correction or 
#'  correction and transformation.
#' @param compared_to the type of data in df_after either "Correction" or 
#'  "Correction and Transformation"
#' @param before_label Character, label for the "before" group in the legend.
#' @param after_label Character, label for the "after" group in the legend.
#'
#' @return A ggplot object with two panels: QC on the left, NonQC on the right.
#' @importFrom dplyr mutate bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap labs theme_bw theme
#'   element_text element_rect element_line scale_fill_manual scale_color_manual
#' @keywords internal
#' @noRd
plot_met_rsd_distributions <- function(df_before,
                                   df_after,
                                   compared_to,
                                   before_label = "Before",
                                   after_label  = "After") {
  
  rsd_before <- metabolite_rsd(df_before)
  rsd_after <- metabolite_rsd(df_after)
  # Basic checks
  required_cols <- c("Metabolite", "RSD_QC", "RSD_NonQC")
  if (!all(required_cols %in% names(rsd_before))) {
    stop("`rsd_before` must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!all(required_cols %in% names(rsd_after))) {
    stop("`rsd_after` must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Tag before/after and combine
  rsd_before2 <- dplyr::mutate(rsd_before, dataset = before_label)
  rsd_after2  <- dplyr::mutate(rsd_after,  dataset = after_label)
  
  rsd_all <- dplyr::bind_rows(rsd_before2, rsd_after2)
  
  # Long format: one column for QC vs NonQC, one for RSD values
  rsd_long <- tidyr::pivot_longer(
    rsd_all,
    cols = c("RSD_QC", "RSD_NonQC"),
    names_to = "type",
    values_to = "RSD"
  )
  
  # QC on left, NonQC on right
  rsd_long$type <- factor(
    rsd_long$type,
    levels = c("RSD_NonQC", "RSD_QC"),
    labels = c("Samples", "QC")
  )
  
  # Drop any NA RSD values
  rsd_long <- rsd_long[!is.na(rsd_long$RSD), , drop = FALSE]
  
  col_vals <- stats::setNames(
    c("#1F77B4", "#FF7F0E"),
    c(before_label, after_label)
  )
  
  p <- ggplot2::ggplot(rsd_long, ggplot2::aes(x = RSD, fill = dataset, color = dataset)) +
    ggplot2::geom_density(alpha = 0.3, adjust = 1) +
    ggplot2::facet_wrap(~ type, nrow = 1, scales = "fixed") +
    ggplot2::labs(
      title = paste("Comparison of RSD Before and After", compared_to),
      x = "RSD (%)",
      y = "Density",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "white", colour = "grey30"),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      plot.title   = ggplot2::element_text(
        size = 14,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title   = ggplot2::element_text(size = 14, face = "bold"),
      axis.text    = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    ) +
    ggplot2::scale_fill_manual(values = col_vals) +
    ggplot2::scale_color_manual(values = col_vals)
  
  p
}

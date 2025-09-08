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
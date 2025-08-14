# rsd comparison plots
source("R/processing_helpers.R")

lab_levels   <- c("Increased","No Change","Decreased")
color_values <- c("Increased"="#B22222",
                  "No Change"="gray25",
                  "Decreased"="#234F1E")

# facet strip labels with percents helper
facet_label_map <- function(df) {
  by_type <- df %>%
    dplyr::group_split(Type) %>%
    purrr::map(~{
      pt <- pct_tbl(.x)
      paste0(
        unique(.x$Type), " — ",
        "Increased ",  pt$percent[pt$change=="Increased"],  "% · ",
        "No change ",  pt$percent[pt$change=="No Change"],   "% · ",
        "Decreased ",  pt$percent[pt$change=="Decreased"],   "%"
      )
    }) %>% unlist()
  stats::setNames(by_type, unique(df$Type))
}

# comparison scatter plot helper
mk_plot <- function(d_all, x, y, facet_labs) {
  if (!nrow(d_all)) return(ggplot2::ggplot() + ggplot2::labs(title = paste0(labs_title, " (no points)")))
  ggplot2::ggplot(d_all, ggplot2::aes(x=.data[[x]], y=.data[[y]], color = change)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = lab_levels,
      labels = c("Increased","No change","Decreased"),
      name   = "RSD Change"
    ) +
    ggplot2::facet_wrap(~ Type, nrow = 1, labeller = ggplot2::as_labeller(facet_labs)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 12),
      axis.text    = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 10),
      legend.background      = ggplot2::element_rect(fill = "white", linewidth = 0.5, linetype = "solid"),
      panel.border           = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::labs(
      title = "Comparison of RSD Before and After Correction",
      x = "RSD Before",
      y = "RSD After"
    )
}

# categorize changes helper function
tag_changes <- function(d, a_before, a_after) {
  if (!nrow(d)) return(dplyr::mutate(d, change = factor(character(), levels = lab_levels)))
  dplyr::mutate(d,
                change = dplyr::case_when(
                  .data[[a_after]]  > .data[[a_before]] ~ "Increased",
                  .data[[a_after]]  < .data[[a_before]] ~ "Decreased",
                  TRUE                                   ~ "No Change"
                ),
                change = factor(change, levels = lab_levels)
  )
}

# Compute percentages
pct_tbl <- function(d) {
  total <- nrow(d); if (!total) return(setNames(data.frame(change=factor(lab_levels, levels=lab_levels),
                                                           percent=c(0,0,0)), c("change","percent")))
  d %>% dplyr::count(change, .drop = FALSE) %>%
    tidyr::complete(change = factor(lab_levels, levels = lab_levels), fill = list(n = 0)) %>%
    dplyr::mutate(percent = round(100 * n / total, 1)) %>% dplyr::select(change, percent)
}

# compute legend percentages
label_map <- function(perc) setNames(
  paste0(c("Increased: ","No change: ","Decreased: "), perc$percent, "%"),
  lab_levels
)

# RSD comparison scatter plots when computing RSD by metabolite
plot_rsd_comparison <- function(df_before, df_after) {

  rsdBefore <- metabolite_rsd(df_before)
  rsdAfter <- metabolite_rsd(df_after)
  
  # Merge data on Metabolite
  df <- dplyr::inner_join(
    dplyr::rename(rsdBefore, rsd_qc_before = RSD_QC, rsd_nonqc_before = RSD_NonQC),
    dplyr::rename(rsdAfter,  rsd_qc_after  = RSD_QC, rsd_nonqc_after  = RSD_NonQC),
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
    dplyr::transmute(Type = "Samples", before = rsd_nonqc_before, after = rsd_nonqc_after, change)
  
  df_qcs <- df %>%
    dplyr::filter(is.finite(rsd_qc_before), is.finite(rsd_qc_after)) %>%
    tag_changes("rsd_qc_before", "rsd_qc_after") %>%
    dplyr::transmute(Type = "QCs", before = rsd_qc_before, after = rsd_qc_after, change)
  
  d_all <- dplyr::bind_rows(df_samples, df_qcs) %>%
    dplyr::mutate(Type = factor(Type, levels = c("Samples","QCs")),
                  change = factor(change, levels = lab_levels))
  
  # make legend map
  facet_labs <- facet_label_map(d_all)
  
  # single faceted ggplot
  p <- ggplot2::ggplot(d_all, ggplot2::aes(x = before, y = after, color = change)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = lab_levels,
      labels = c("Increased","No change","Decreased"),
      name   = "RSD Change"
    ) +
    ggplot2::facet_wrap(~ Type, nrow = 1, labeller = ggplot2::as_labeller(facet_labs)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title   = ggplot2::element_text(size = 12),
      axis.text    = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 10),
      legend.background      = ggplot2::element_rect(fill = "white", linewidth = 0.5, linetype = "solid"),
      panel.border           = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::labs(
      title = "Comparison of RSD Before and After Correction",
      x = "RSD Before",
      y = "RSD After"
    )
  
  return(p)
}

# RSD comparison scatter plot when computing RSD by metabolite-class
plot_rsd_comparison_class_met <- function(df_before, df_after) {
  rsdBefore <- class_metabolite_rsd(df_before)
  rsdAfter <- class_metabolite_rsd(df_after)
  
  # Merge data on class and Metabolite
  df <- dplyr::inner_join(
    dplyr::rename(rsdBefore, rsd_before = RSD),
    dplyr::rename(rsdAfter,  rsd_after  = RSD),
    by = c("class", "Metabolite")
  )
  # If empty, return empty plot
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No overlapping metabolites"))
  }
  
  df_samples <- df[df$class != "QC", ]
  df_qcs <- df[df$class == "QC", ]
  
  # Get only finite data and Categorize changes
  df_samples <- df_samples %>%
    dplyr::filter(is.finite(rsd_before), is.finite(rsd_after)) %>%
    tag_changes("rsd_before", "rsd_after") %>%
    dplyr::transmute(Type = "Samples", before = rsd_before, after = rsd_after, change)
  
  df_qcs <- df_qcs %>%
    dplyr::filter(is.finite(rsd_before), is.finite(rsd_after)) %>%
    tag_changes("rsd_before", "rsd_after") %>%
    dplyr::transmute(Type = "QCs", before = rsd_before, after = rsd_after, change)
  
  d_all <- dplyr::bind_rows(df_samples, df_qcs) %>%
    dplyr::mutate(Type = factor(Type, levels = c("Samples","QCs")),
                  change = factor(change, levels = lab_levels))
  
  # make legend map
  facet_labs <- facet_label_map(d_all)
  
  # single faceted ggplot
  mk_plot(d_all, "before", "after", facet_labs)
}
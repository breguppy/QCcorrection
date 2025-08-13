# rsd comparison plots

#library(ggplot2)
#library(patchwork)
#library(dplyr)
#library(tidyr)
#library(tidyverse)
#library(grid)
#library(gridExtra)
source("R/processing_helpers.R")

lab_levels   <- c("Increased","No Change","Decreased")
color_values <- c("Increased"="#B22222",
                  "No Change"="gray25",
                  "Decreased"="#234F1E")


# comparison scatter plot helper
mk_plot <- function(d, x, y, labs_title) {
  if (!nrow(d)) return(ggplot2::ggplot() + ggplot2::labs(title = paste0(labs_title, " (no points)")))
  ggplot2::ggplot(d, ggplot2::aes(x=.data[[x]], y=.data[[y]], color=change)) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed") +
    ggplot2::geom_point(size=3, na.rm=TRUE) +
    ggplot2::scale_color_manual(values=color_values, breaks=names(labs_title), labels=labs_title, name="RSD Change") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size=14, hjust=0.5, face="bold"),
      axis.title = ggplot2::element_text(size=12),
      axis.text  = ggplot2::element_text(size=10),
      legend.text = ggplot2::element_text(size=10),
      legend.position = "inside",
      legend.position.inside = c(0, 1),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill="white", linewidth=0.5, linetype="solid"),
      panel.border = ggplot2::element_rect(color="black", fill=NA, linewidth=1)
    ) +
    ggplot2::labs(x="RSD Before", y="RSD After")
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
  
  print("IN PLOTTING FUNCTION")
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
  
  print("DATA MERGED")
  # Get only finite data:
  df_samples <- dplyr::filter(df, is.finite(rsd_nonqc_before), is.finite(rsd_nonqc_after))
  df_qcs     <- dplyr::filter(df, is.finite(rsd_qc_before),    is.finite(rsd_qc_after))
  
  # categorize changes
  df_samples <- tag_changes(df_samples, "rsd_nonqc_before", "rsd_nonqc_after")
  df_qcs     <- tag_changes(df_qcs,     "rsd_qc_before",    "rsd_qc_after")
  
  # make legend map
  label_samples <- label_map(pct_tbl(df_samples))
  label_qcs     <- label_map(pct_tbl(df_qcs))
  
  print("DATA CATEGORIZED")
  # before and after plots
  p1 <- mk_plot(df_samples, "rsd_nonqc_before", "rsd_nonqc_after", label_samples) + ggplot2::labs(title="Samples")
  p2 <- mk_plot(df_qcs,     "rsd_qc_before",    "rsd_qc_after",    label_qcs)     + ggplot2::labs(title="QCs")
  
  print("INDIVIDUAL PLOTS MADE")
  # patch them together
  plot <- patchwork::wrap_plots(p1, p2, nrow = 1) +
    patchwork::plot_annotation(
      "Comparison of RSD Before and After Correction",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  print("COMBINED PLOT MADE")
  
  return(plot)
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
  
  # Get only finite data:
  df_samples <- dplyr::filter(df_samples, is.finite(rsd_before), is.finite(rsd_after))
  df_qcs     <- dplyr::filter(df_qcs, is.finite(rsd_before),    is.finite(rsd_after))
  
  # Categorize changes in CV
  df_samples <- tag_changes(df_samples, "rsd_before", "rsd_after")
  df_qcs     <- tag_changes(df_qcs,     "rsd_before",    "rsd_after")
  
  # make legend map
  label_samples <- label_map(pct_tbl(df_samples))
  label_qcs     <- label_map(pct_tbl(df_qcs))
  
  # before and after plots
  p1 <- mk_plot(df_samples, "rsd_before", "rsd_after", label_samples) + ggplot2::labs(title="Samples")
  p2 <- mk_plot(df_qcs,     "rsd_before",    "rsd_after",    label_qcs)     + ggplot2::labs(title="QCs")
  
  # patch them together
  patchwork::wrap_plots(p1, p2, nrow = 1) +
    patchwork::plot_annotation(
      "Comparison of RSD Before and After Correction",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5))
    )
}
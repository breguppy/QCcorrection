# Matabolite scatter plot for random forest corrected data

met_scatter_rf <- function(data_raw, data_cor, i) {
  # add panel tags
  data_raw <- dplyr::mutate(data_raw, type = ifelse(class=="QC","QC","Sample"),
                            panel = factor("Raw", levels=c("Raw","Corrected")))
  data_cor <- dplyr::mutate(data_cor, type = ifelse(class=="QC","QC","Sample"),
                            panel = factor("Corrected", levels=c("Raw","Corrected")))
  df_all   <- dplyr::bind_rows(data_raw, data_cor)
  
  # qc stats
  get_stats <- function(df, panel) {
    qc_vals <- df[[i]][df$type=="QC"]
    m <- mean(qc_vals, na.rm=TRUE)
    s <- sd(qc_vals, na.rm=TRUE)
    tibble::tibble(
      y=c(m+s,m-s,m+2*s,m-2*s),
      sd=factor(c("\u00B1 1 SD","\u00B1 1 SD","\u00B1 2 SD","\u00B1 2 SD"), levels=c("\u00B1 1 SD","\u00B1 2 SD")),
      panel=factor(panel, levels=c("Raw","Corrected"))
    )
  }
  sd_df <- dplyr::bind_rows(get_stats(data_raw,"Raw"),
                            get_stats(data_cor,"Corrected"))
  
  # batch shading
  get_batches <- function(df,panel) {
    df |>
      dplyr::group_by(batch) |>
      dplyr::summarize(xmin=min(order), xmax=max(order), .groups="drop") |>
      dplyr::arrange(xmin) |>
      dplyr::mutate(fill=rep(c("lightgray","white"), length.out=dplyr::n()),
                    panel=factor(panel, levels=c("Raw","Corrected")))
  }
  batch_ranges <- dplyr::bind_rows(get_batches(data_raw,"Raw"),
                                   get_batches(data_cor,"Corrected"))
  
  color_scale <- ggplot2::scale_color_manual(name="Type:",
                                             values=c(Sample="#F5C710", QC="#305CDE"))
  lty_scale <- ggplot2::scale_linetype_manual(name="SD Range:",
                                              values=c("\u00B1 1 SD"="dashed","\u00B1 2 SD"="solid"),
                                              guide=ggplot2::guide_legend(
                                                override.aes=list(color=c("grey20","#950606"))
                                              ))
  
  ggplot2::ggplot(df_all, ggplot2::aes(x=order, y=.data[[i]])) +
    ggplot2::geom_rect(data=batch_ranges,
                       ggplot2::aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=fill),
                       inherit.aes=FALSE, alpha=0.3, show.legend=FALSE) +
    ggplot2::scale_fill_identity(guide="none") +
    ggplot2::geom_hline(data=sd_df |> dplyr::filter(sd=="\u00B1 1 SD"),
                        ggplot2::aes(yintercept=y, linetype=sd),
                        color="grey20", linewidth=0.75) +
    ggplot2::geom_hline(data=sd_df |> dplyr::filter(sd=="\u00B1 2 SD"),
                        ggplot2::aes(yintercept=y, linetype=sd),
                        color="#950606", linewidth=0.75) +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "Sample"),
      ggplot2::aes(color = type), size = 2
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_all, type == "QC"),
      ggplot2::aes(color = type), size = 2
    ) +
    color_scale + lty_scale +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(title=i, x="Injection Order", y="Intensity") +
    ggplot2::theme_minimal(base_size=10) +
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
      legend.key.width  = grid::unit(0.5, "cm"),   # shrink keys
      legend.key.height = grid::unit(0.3, "cm"),
      legend.margin     = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      legend.box.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2),
      strip.text.x = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.5, 0.5, 0.8, 0.5), "cm")
    )
}


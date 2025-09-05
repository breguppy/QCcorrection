#' exports figures and returns file path for zip folder
#'
#' @keywords internal
#' @noRd
export_figures <- function(p, d, out_dir = tempdir()) {
  .require_pkg("ggplot2", "write figures")
  fig_dir <- file.path(out_dir, "figures")
  if (dir.exists(fig_dir))
    unlink(fig_dir, recursive = TRUE)
  dir.create(fig_dir)
  fmt <- match.arg(p$fig_format, c("png", "pdf"))
  
  mk <- function(subdir) {
    x <- file.path(fig_dir, subdir)
    if (!dir.exists(x))
      dir.create(x, recursive = TRUE)
    x
  }
  rsd_dir <- mk("RSD figures")
  pca_dir <- mk("PCA plots")
  met_dir <- mk("metabolite figures")
  
  save_plot <- function(path, plot, w, h) {
    if (fmt == "png") {
      ggplot2::ggsave(
        path,
        plot = plot,
        width = w,
        height = h,
        units = "in",
        dpi = 300,
        bg = "white"
      )
    } else {
      ggplot2::ggsave(
        path,
        plot = plot,
        width = w,
        height = h,
        units = "in",
        device = grDevices::cairo_pdf
      )
    }
    normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  
  raw_cols <- setdiff(names(d$filtered$df), c("sample", "batch", "class", "order"))
  cor_cols <- setdiff(names(d$filtered_corrected$df),
                      c("sample", "batch", "class", "order"))
  cols <- intersect(raw_cols, cor_cols)
  met_paths <- character(0)
  n <- length(cols)
  N <- n + 2
  shiny::withProgress(message = "Creating figures...", value = 0, {
    rsd_plot <- make_rsd_plot(p, d)
    rsd_path <- file.path(rsd_dir, sprintf("rsd_comparison_%s.%s", p$rsd_cal, fmt))
    rsd_path <- save_plot(rsd_path, rsd_plot, 7.5, 4.5)
    shiny::incProgress(1 / N, detail = "Saved: rsd figure")
    
    pca_plot <- make_pca_plot(p, d)  # ensure this uses p$color_col safely
    pca_path <- file.path(pca_dir,
                          sprintf("pca_comparison_%s.%s", p$color_col, fmt))
    pca_path <- save_plot(pca_path, pca_plot, 8.333, 4.417)
    shiny::incProgress(1 / N, detail = "Saved: pca figure")
    
    for (i in seq_len(n)) {
      metab <- cols[i]
      fig <- make_met_scatter(d, metab)
      safe <- gsub("[^A-Za-z0-9_\\-]+", "_", metab)
      path <- file.path(met_dir, sprintf("%s.%s", safe, fmt))
      met_paths <- c(met_paths, save_plot(path, fig, 5, 5))
      shiny::incProgress(1 / n, detail = paste("Saved:", safe))
    }
  })
  
  list(
    fig_dir = normalizePath(fig_dir, winslash = "/", mustWork = TRUE),
    rsd = rsd_path,
    pca = pca_path,
    metabolite = met_paths
  )
}

#' Zips figure folder
#'
#' @keywords internal
#' @noRd
zip_figures_dir <- function(fig_dir, zip_file) {
  .require_pkg("zip", "create a zip archive")
  fig_dir <- normalizePath(fig_dir, winslash = "/", mustWork = TRUE)
  zip::zipr(zip_file, files = fig_dir)   
  zip_file
}

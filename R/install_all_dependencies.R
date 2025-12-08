#' Install all optional dependencies
#'
#' Installs all packages listed in the \strong{Suggests} field of QCcorrection that
#' enable optional functionality such as random forest correction, outlier detection,
#' reporting, testing, and PDF export.
#'
#' @details
#' The function checks which of the optional packages are not already installed
#' and installs only the missing ones from CRAN. It does not reinstall packages
#' that are already present.
#'
#' This helper is useful for users installing QCcorrection from GitHub with
#' \code{remotes::install_github()}, which does not automatically install
#' packages listed in \emph{Suggests}.
#'
#' The optional packages are:
#' \itemize{
#'   \item \pkg{randomForest}, \pkg{robustbase}, \pkg{outliers}, \pkg{EnvStats}, \pkg{corpcor} – statistical methods for correction and outlier detection
#'   \item \pkg{ggtext}, \pkg{cowplot} – enhanced plotting and figure assembly
#'   \item \pkg{pagedown}, \pkg{rmarkdown}, \pkg{zip}, \pkg{knitr} – report generation and export
#'   \item \pkg{httpuv}, \pkg{jsonlite}, \pkg{pkgload} – internal testing and Shiny app support
#'   \item \pkg{shinytest2}, \pkg{testthat}, \pkg{chromote} – automated testing and headless browser rendering
#' }
#'
#' @return Invisibly returns a character vector of newly installed package names.
#' @examples
#' \dontrun{
#' QCcorrection::install_all_dependencies()
#' }
#' @export
install_all_dependencies <- function() {
  
  # ---- CRAN optional packages ----
  cran_pkgs <- c(
    "randomForest", "robustbase", "outliers", "EnvStats",
    "ggtext", "cowplot", "pagedown", "rmarkdown", "zip",
    "corpcor", "httpuv", "jsonlite", "shinytest2",
    "testthat", "chromote", "pkgload", "knitr"
  )
  
  # Identify missing CRAN packages
  cran_missing <- setdiff(cran_pkgs, rownames(installed.packages()))
  
  
  # ---- Bioconductor packages ----
  bioc_pkgs <- c("impute")
  
  # Install BiocManager if needed
  if (!("BiocManager" %in% rownames(installed.packages()))) {
    install.packages("BiocManager")
    message("Installed BiocManager")
  }
  
  # Identify missing Bioconductor packages
  bioc_missing <- setdiff(
    bioc_pkgs,
    rownames(installed.packages())
  )
  
  # Install missing Bioconductor packages
  if (length(bioc_missing) > 0) {
    BiocManager::install(bioc_missing, ask = FALSE)
    message("Installed Bioconductor packages: ",
            paste(bioc_missing, collapse = ", "))
  }
  
  
  # ---- Install missing CRAN packages ----
  if (length(cran_missing) > 0) {
    install.packages(cran_missing)
    message("Installed CRAN packages: ",
            paste(cran_missing, collapse = ", "))
  }
  
  # ---- Final message ----
  if (length(cran_missing) == 0 && length(bioc_missing) == 0) {
    message("All optional QCcorrection dependencies already installed.")
  }
  
  invisible(list(
    cran = cran_missing,
    bioc = bioc_missing
  ))
}
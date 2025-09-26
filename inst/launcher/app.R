# app.R for creating desktop app
#' @keywords internal
#' @noRd
#' 
paths <- c(
"C:/Program Files/Quarto/bin/tools",
"C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
"C:/Program Files/RStudio/bin/pandoc"
)
p <- paths[file.exists(file.path(paths, "pandoc.exe"))][1]
if (!is.na(p)) Sys.setenv(RSTUDIO_PANDOC = p)

# quick assert for visibility during troubleshooting
if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available()) {
  stop("Pandoc not available. RSTUDIO_PANDOC=", Sys.getenv("RSTUDIO_PANDOC"))
}
if (!requireNamespace("QCcorrection", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  
  # Absolute path to your package root. Note the slash after C:\.
  dev_path <- "C:/Users/bguppy/Documents/QCcorrection"
  dev_path <- normalizePath(dev_path, winslash = "/", mustWork = TRUE)
  
  remotes::install_local(dev_path, upgrade = "never", force = TRUE)
}
if (Sys.getenv("RSTUDIO_PANDOC") == "" && file.exists("C:/Program Files/Quarto/bin/tools/pandoc.exe")) {
  Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/Quarto/bin/tools")
}

QCcorrection::run_app()
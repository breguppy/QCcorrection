#!/usr/bin/env Rscript
# Use only the bundled library
lib_dir <- Sys.getenv("R_PACK_LIB")
if (nzchar(lib_dir)) .libPaths(c(normalizePath(lib_dir), .libPaths()))

# Optional: bundled pandoc
pandoc <- Sys.getenv("PANDOC_PATH")
if (nzchar(pandoc)) Sys.setenv(RSTUDIO_PANDOC = normalizePath(pandoc))

# Host/port from Electron
port <- as.integer(Sys.getenv("SHINY_PORT", "3170"))
options(shiny.host = "127.0.0.1", shiny.port = port, shiny.launch.browser = FALSE)

# Load your package and run
library(QCcorrection)
run_app()
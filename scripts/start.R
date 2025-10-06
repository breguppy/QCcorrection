#!/usr/bin/env Rscript

# 1) library path
lib_dir <- Sys.getenv("R_PACK_LIB", "C:/Users/bguppy/Documents/QCcorrection/r-env/win/library")
.libPaths(c(normalizePath(lib_dir), .libPaths()))

# 2) optional Pandoc
pandoc_dir <- "C:/Users/bguppy/Documents/QCcorrection/r-env/win/pandoc/pandoc-3.8.2"
if (dir.exists(pandoc_dir)) {
  Sys.setenv(RSTUDIO_PANDOC = normalizePath(pandoc_dir))
}

# 3) Shiny host/port
port <- as.integer(Sys.getenv("SHINY_PORT", "3170"))
options(shiny.host = "127.0.0.1", shiny.port = port, shiny.launch.browser = FALSE)

# 4) run app
library(QCcorrection)
run_app()

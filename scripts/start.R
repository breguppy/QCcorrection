#!/usr/bin/env Rscript

# 1) library path
lib_dir <- Sys.getenv("R_PACK_LIB")
if (nzchar(lib_dir)) .libPaths(c(normalizePath(lib_dir), .libPaths()))

# 2) optional Pandoc
pd <- Sys.getenv("RSTUDIO_PANDOC")
if (!nzchar(pd)) {
  cand <- c(
    "r-env/win/pandoc",
    "r-env/win/pandoc/pandoc-3.8.2"
  )
  hit <- cand[file.exists(file.path(cand, "pandoc.exe"))][1]
  if (!is.na(hit)) Sys.setenv(RSTUDIO_PANDOC = normalizePath(hit))
}

# 3) Shiny host/port
port <- as.integer(Sys.getenv("SHINY_PORT", "3170"))
options(shiny.host = "127.0.0.1", shiny.port = port, shiny.launch.browser = FALSE)
Sys.setenv(SHINY_SERVER_VERSION = "0.3.4")

# 4) run app
library(QCcorrection)
app <- run_app()
shiny::runApp(app, host = "127.0.0.1", port = port, launch.browser = FALSE)

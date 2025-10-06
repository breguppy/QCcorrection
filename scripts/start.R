#!/usr/bin/env Rscript

# 1) library path
lib_dir <- Sys.getenv("R_PACK_LIB")
if (nzchar(lib_dir)) .libPaths(c(normalizePath(lib_dir), .libPaths()))

# 2) optional Pandoc
pd <- Sys.getenv("RSTUDIO_PANDOC"); if (nzchar(pd)) Sys.setenv(RSTUDIO_PANDOC = normalizePath(pd))

# 3) Shiny host/port
port <- as.integer(Sys.getenv("SHINY_PORT", "3170"))
options(shiny.host = "127.0.0.1", shiny.port = port, shiny.launch.browser = FALSE)

# 4) run app
library(QCcorrection)
run_app()

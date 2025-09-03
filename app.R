#' Launches QCcorrection Shiny App
#' @export
# source("R/app_ui.R"); source("R/app_server.R")
run_app <- function() shinyApp(ui = app_ui(), server = app_server)
run_app()

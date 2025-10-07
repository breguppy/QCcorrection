#' Launches QCcorrection Shiny App
#' @export

run_app <- function() {
  options(shiny.launch.browser = FALSE, shiny.testmode = FALSE)
  if (Sys.getenv("ELECTRON", "0") != "1") {
    options(shiny.launch.browser = TRUE)
  }
  shiny::shinyApp(ui = app_ui(), server = app_server)
}

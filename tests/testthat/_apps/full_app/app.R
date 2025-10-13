options(shiny.testmode = TRUE)
library(shiny)
if (!requireNamespace("QCcorrection", quietly = TRUE)) {
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(QCcorrection)

shinyApp(
  ui     = QCcorrection:::app_ui(),
  server = function(input, output, session) QCcorrection:::app_server(input, output, session)
)
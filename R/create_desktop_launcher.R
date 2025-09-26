#' Create a desktop launcher directory and .bat via shiny.exe
#' @export
create_desktop_launcher <- function(path = "QCcorrection_app", app_name = "QCcorrection") {
  if (!requireNamespace("shiny.exe", quietly = TRUE)) {
    stop("Install the 'shiny.exe' package first.")
  }
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  src <- system.file("launcher", "app.R", package = "QCcorrection")
  if (nzchar(src)) file.copy(src, file.path(path, "app.R"), overwrite = TRUE)
  shiny.exe::shiny.exe(appName = app_name, appDir = normalizePath(path))
  invisible(normalizePath(path))
}

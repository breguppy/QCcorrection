#' Create a desktop launcher directory and .bat via shiny.exe
#' @export

make_launcher_bat <- function(dir = system.file("launcher", package = "QCcorrection"),
                              bat_name = "QCcorrection.bat",
                              app_file = "app.R") {
  rs <- normalizePath(file.path(R.home("bin"), "Rscript.exe"), winslash = "\\")
  bat <- sprintf(
    '@echo off
setlocal
set Rscript="%s"
set APPDIR=%%~dp0
%%Rscript%% --vanilla "%%APPDIR%%%s"
if errorlevel 1 (
  echo Failed to start QCcorrection app. Press any key to exit.
  pause >nul
)
endlocal',
  rs, paste0("\\", app_file)
  )
writeLines(bat, file.path(dir, bat_name))
}
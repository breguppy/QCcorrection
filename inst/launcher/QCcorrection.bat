@echo off
setlocal
set Rscript="C:\Users\bguppy\AppData\Local\Programs\R\R-4.4.1\bin\x64\Rscript.exe"
set APPDIR=%~dp0
%Rscript% --vanilla "%APPDIR%\app.R"
if errorlevel 1 (
  echo Failed to start QCcorrection app. Press any key to exit.
  pause >nul
)
endlocal

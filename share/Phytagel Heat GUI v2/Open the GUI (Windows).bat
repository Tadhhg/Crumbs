@echo off
rem  Opens the Phytagel heat GUI.  Nothing to install, nothing to type.
rem  Needs MATLAB R2022b or newer already on this machine.
setlocal
set "CODE=%~dp0code\Heatmodel\matlab"

if not exist "%CODE%\heat_gui.m" (
  echo.
  echo   Cannot find the model code next to this file.
  echo   Keep this .bat inside the "Phytagel Heat GUI" folder -- it needs
  echo   the "code" folder that sits beside it.
  echo.
  pause & exit /b 1
)

echo Starting MATLAB and opening the GUI -- this takes about 30 seconds...

where matlab >nul 2>&1
if %errorlevel%==0 (
  start "" matlab -sd "%CODE%" -r "heat_gui"
  exit /b 0
)

rem MATLAB is not on the PATH -- look for the newest install instead.
for %%P in ("%ProgramFiles%\MATLAB" "%ProgramFiles(x86)%\MATLAB" "C:\MATLAB") do (
  if exist %%P (
    for /f "delims=" %%D in ('dir /b /o-n %%P 2^>nul') do (
      if exist "%%~P\%%D\bin\matlab.exe" (
        start "" "%%~P\%%D\bin\matlab.exe" -sd "%CODE%" -r "heat_gui"
        exit /b 0
      )
    )
  )
)

echo.
echo   Could not find MATLAB on this computer.
echo   Install MATLAB first (your university almost certainly provides it),
echo   then double-click this file again.
echo.
pause

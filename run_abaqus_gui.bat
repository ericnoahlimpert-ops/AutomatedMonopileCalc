@echo off
rem Opens Abaqus/CAE and runs the automated monopile model builder (main.py)
rem there; the model stays open afterwards. This window shows the progress
rem (finished load levels and calculations) live. Without GUI:
rem run_abaqus.bat. The progress is also saved in RESULTS_DIR\progress.log.
rem
rem Usage:
rem   run_abaqus_gui.bat                  input folder from INPUT_DIR below
rem   run_abaqus_gui.bat <input folder>   e.g. Inputs\Verification\PISA_Verification\C1_ultimate
rem   (or drag an input folder onto this file)

setlocal

rem Settings
rem   INPUT_DIR        folder with Input.py, AdvancedSettings.py and the input
rem                    CSV files (relative to this folder or absolute);
rem                    a folder given as argument overrides it
rem   STAY_OPEN=TRUE   keep this window open at the end until a key is pressed
rem   STAY_OPEN=FALSE  close it at once (e.g. when called from another script)
set INPUT_DIR=Inputs
set STAY_OPEN=TRUE

if not "%~1"=="" set "INPUT_DIR=%~1"
cd /d "%~dp0"

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_abaqus.ps1" -Mode script -InputDir "%INPUT_DIR%"
set EXITCODE=%ERRORLEVEL%

if %EXITCODE% neq 0 (
    echo.
    echo Abaqus/CAE closed with an error ^(exit code %EXITCODE%^).
)

rem Keep the window open (STAY_OPEN, see top)
if /i "%STAY_OPEN%"=="TRUE" pause
endlocal & exit /b %EXITCODE%

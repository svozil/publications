@echo off
setlocal

:: ============================================================================
::  get.bat (Get a specific file from the GitHub repository)
::
::  This script simplifies downloading a single file. You only need to
::  provide the filename.
::
::  Wildcards (*) are not supported for downloads. If you want to get
::  all new or updated files, you should use the "update.bat" script
::  (which runs "git pull").
:: ============================================================================

:: --- Configuration ---
set "GITHUB_USER=svozil"
set "GITHUB_REPO=publications"
set "BRANCH=master"
:: ---------------------

:: Construct the base URL from the configuration above
set "BASE_URL=https://raw.githubusercontent.com/%GITHUB_USER%/%GITHUB_REPO%/%BRANCH%"

:: --- Script Logic ---

:: Check if a filename was provided
if "%~1"=="" (
    echo [ERROR] You must provide a filename.
    echo.
    echo   Usage for a single file:
    echo     get "path/to/your-file.pdf"
    echo.
    echo   To get ALL new files, run the update.bat script instead.
    goto :eof
)

:: Check if the user tried to use a wildcard
echo "%~1" | findstr "*" >nul
if errorlevel 1 (
    goto :DownloadSingleFile
) else (
    goto :HandleWildcard
)

:DownloadSingleFile
    set "FILENAME=%~1"
    set "FULL_URL=%BASE_URL%/%FILENAME%"

    echo.
    echo [INFO] Downloading single file:
    echo   %FULL_URL%
    echo.
    echo [INFO] Saving as:
    echo   %FILENAME%
    echo -------------------------------------------------

    :: Use curl to download the specific file
    curl -L -o "%FILENAME%" "%FULL_URL%"

    echo -------------------------------------------------
    echo.
    echo [SUCCESS] Download complete.
    goto :eof

:HandleWildcard
    echo.
    echo [INFO] Wildcard (*) detected.
    echo.
    echo This script cannot download multiple files with a wildcard.
    echo The correct way to get all new and updated files from GitHub
    echo is to "pull" all the latest changes.
    echo.
    echo To do this, please run the 'update.bat' script.
    echo.
    goto :eof

endlocal
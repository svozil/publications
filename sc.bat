@echo off
setlocal

:: ============================================================================
::  sc.bat (Sync Changes)
::
::  This batch file stages all changes, commits them with a provided
::  message, and pushes them to the remote repository.
::
::  USAGE:
::  From your repository's command line, type:
::  sc "Your descriptive commit message here"
:: ============================================================================

:: Check if a commit message argument was provided.
::if "%~1"=="" (
::    echo.
::    echo [ERROR] No commit message provided.
::    echo.
::    echo   You must provide a reason for the change in quotes.
::    echo   Example: sc "Updated the introduction for my latest paper"
::    echo.
::    goto :eof
::)

echo.
echo [STEP 1/3] Staging all new and modified files...
git add .
echo -------------------------------------------------

echo.
echo [STEP 2/3] Committing changes with your message...
git commit -m "saving new file versions and files"

:: The commit command fails if there are no staged changes.
:: We check for that error, but continue anyway, in case a previous
:: commit failed to push.
if errorlevel 1 (
    echo   (Note: This is normal if there were no new changes to commit.)
)
echo -------------------------------------------------

echo.
echo [STEP 3/3] Pushing changes to GitHub...
git push
echo -------------------------------------------------

echo.
echo Script finished.
endlocal

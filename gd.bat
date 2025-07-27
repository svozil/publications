@echo off
setlocal

:: ============================================================================
::  update.bat (Download from Remote)
::
::  This batch file downloads (pulls) the latest changes from the
::  remote GitHub repository into your local copy.
::
::  Run this before you start making your own changes to ensure you
::  are working with the most up-to-date version.
:: ============================================================================

echo.
echo [INFO] Attempting to download latest changes from GitHub...
echo -------------------------------------------------

:: The "git pull" command fetches and merges changes from the remote.
git pull

echo -------------------------------------------------
echo.
echo [SUCCESS] Your local repository is now up-to-date.
echo.

endlocal
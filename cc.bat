@echo off

:: 1. Copy all matching files in a single pass (separate filters by spaces)
robocopy "C:\mytex" "C:\github\publications" "20*.bbl" "20*.pdf" "20*.tex" "20*.zip" "svozil.bib" /XO

:: 2. Navigate to the Git repository
cd /d "C:\github\publications"

:: 3. Run Git commands (Assuming you meant to stage, commit, and push)
git add .
git commit -m "Auto-update publications"
git push

:: 4. Return to the original directory
cd /d "C:\mytex"

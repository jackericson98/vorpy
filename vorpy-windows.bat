@echo off
setlocal
cd /d "%~dp0"

if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" -m vorpy.workbench %*
) else (
    where py >nul 2>nul
    if not errorlevel 1 (
        py -3 -m vorpy.workbench %*
    ) else (
        python -m vorpy.workbench %*
    )
)

if errorlevel 1 (
    echo.
    echo VorPy could not start. From this directory, install it with:
    echo   py -3 -m pip install -e ".[gui]"
    echo.
    pause
    exit /b 1
)

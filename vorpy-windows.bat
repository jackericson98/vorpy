@echo off
setlocal
cd /d "%~dp0"

set "VORPY_PYTHON=.venv\Scripts\python.exe"

if not exist "%VORPY_PYTHON%" (
    echo ============================================================
    echo VorPy first-time setup
    echo ============================================================
    echo Creating a local Python environment. The dependency download
    echo may take several minutes. This normally happens only once.
    echo.

    where py >nul 2>nul
    if not errorlevel 1 (
        py -3 -m venv .venv
    ) else (
        where python >nul 2>nul
        if errorlevel 1 goto :python_missing
        python -m venv .venv
    )

    if errorlevel 1 goto :setup_failed
)

"%VORPY_PYTHON%" -c "import PySide6, pyvista, pyvistaqt, vorpy.workbench" >nul 2>nul
if errorlevel 1 (
    echo Installing VorPy and its graphical dependencies...
    echo Please keep this window open. This may take several minutes.
    echo.
    "%VORPY_PYTHON%" -m pip install -e ".[gui]"
    if errorlevel 1 goto :setup_failed
    echo.
    echo Installation complete. Starting VorPy...
    echo.
)

"%VORPY_PYTHON%" -m vorpy.workbench %*
if errorlevel 1 goto :launch_failed
exit /b 0

:python_missing
echo.
echo Python 3 was not found. Install Python 3 from https://python.org,
echo then open this launcher again.
echo During installation, enable the option to add Python to PATH.
pause
exit /b 1

:setup_failed
echo.
echo VorPy setup failed. Review the installation error above, then try again.
pause
exit /b 1

:launch_failed
echo.
echo VorPy was installed but could not start. Review the error above.
pause
exit /b 1

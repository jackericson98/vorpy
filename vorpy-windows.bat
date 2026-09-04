@echo off
setlocal
cd /d "%~dp0"

set "VORPY_PYTHON=.venv\Scripts\python.exe"
set "VORPY_PYTHONW=.venv\Scripts\pythonw.exe"
set "VORPY_READY=.venv\.vorpy-gui-ready"

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

if not exist "%VORPY_READY%" (
    "%VORPY_PYTHON%" -c "import importlib.util, sys; names = ('PySide6', 'pyvista', 'pyvistaqt', 'vorpy.workbench'); sys.exit(0 if all(importlib.util.find_spec(name) for name in names) else 1)" >nul 2>nul
)
if not exist "%VORPY_READY%" if errorlevel 1 (
    if exist "%VORPY_PYTHONW%" (
        start "" /wait "%VORPY_PYTHONW%" "vorpy\workbench\bootstrap.py"
    ) else (
        "%VORPY_PYTHON%" "vorpy\workbench\bootstrap.py"
    )
    if errorlevel 1 goto :setup_failed
    type nul > "%VORPY_READY%"
    echo.
    echo Installation complete. Starting VorPy...
    echo.
)

if not exist "%VORPY_READY%" type nul > "%VORPY_READY%"

if exist "%VORPY_PYTHONW%" (
    start "" "%VORPY_PYTHONW%" -m vorpy.workbench %*
) else (
    start "" "%VORPY_PYTHON%" -m vorpy.workbench %*
)
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

@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo === BELLO Setup ===
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python is not installed or not in PATH.
    echo   Download from https://www.python.org/downloads/
    echo   Make sure to check "Add Python to PATH" during installation.
    pause
    exit /b 1
)

python -c "import tkinter" >nul 2>&1
if errorlevel 1 (
    echo [ERROR] tkinter is not available.
    echo   Re-install Python and make sure "tcl/tk and IDLE" is checked.
    pause
    exit /b 1
)
echo [OK] tkinter is available

if not exist ".venv" (
    echo Creating virtual environment...
    python -m venv .venv
)

echo Activating virtual environment...
call .venv\Scripts\activate.bat

echo Installing dependencies...
pip install --upgrade pip -q
pip install -r requirements.txt -q

echo.
echo === Setup complete ===
echo.
echo To run BELLO:
echo   .venv\Scripts\activate.bat
echo   python BELLO_GUI.py
echo.
pause

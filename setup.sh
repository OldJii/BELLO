#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== BELLO Setup ==="
echo ""

# Detect Python version for tkinter package name
PY_VER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

# Auto-install tkinter if missing
if ! python3 -c "import tkinter" 2>/dev/null; then
    echo "[INFO] tkinter not found. Attempting to install..."
    if command -v brew &>/dev/null; then
        echo "  -> brew install python-tk@${PY_VER}"
        brew install "python-tk@${PY_VER}" || true
    elif command -v apt &>/dev/null; then
        echo "  -> sudo apt install python3-tk"
        sudo apt install -y python3-tk || true
    elif command -v dnf &>/dev/null; then
        echo "  -> sudo dnf install python3-tkinter"
        sudo dnf install -y python3-tkinter || true
    else
        echo "[WARNING] Cannot auto-install tkinter."
        echo "  Please install it manually for your platform."
    fi
fi

# Verify tkinter
if python3 -c "import tkinter" 2>/dev/null; then
    echo "[OK] tkinter is available"
else
    echo "[ERROR] tkinter is still not available. GUI will not work."
    echo "  macOS:   brew install python-tk@${PY_VER}"
    echo "  Ubuntu:  sudo apt install python3-tk"
    echo "  Fedora:  sudo dnf install python3-tkinter"
    exit 1
fi

# Create / rebuild venv
if [ ! -d ".venv" ] || [ ! -f ".venv/bin/python" ]; then
    echo "Creating virtual environment..."
    rm -rf .venv
    python3 -m venv .venv
fi

echo "Activating virtual environment..."
source .venv/bin/activate

echo "Installing Python dependencies..."
pip install --upgrade pip -q
pip install -r requirements.txt -q

echo ""
echo "=== Setup complete ==="
echo ""
echo "To run BELLO:"
echo "  source .venv/bin/activate"
echo "  python BELLO_GUI.py"

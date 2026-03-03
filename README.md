# BELLO

BELLO (Bond Element Lattice Locality Order) is a post-processing toolkit for characterizing local structural order in amorphous and disordered materials. It works with standard XYZ trajectory files produced by any atomistic simulation code (DFT, MD, etc.) and provides automated identification of coordination environments, order-parameter statistics, radial/angular distribution analysis, and publication-ready visualizations — all through a cross-platform graphical interface with Chinese/English language support.

## Features

- **Local-order classification** — automatically identifies 0- through 6-fold coordination, tetrahedral, and octahedral motifs for every atom in every frame
- **Order parameter *q*** — computes the Steinhardt-type bond-orientational order parameter to quantify structural regularity
- **Automatic threshold** — determines the nearest-neighbor cutoff from the first peak of the radial distribution function via Gaussian fitting
- **Radial pair distribution function** — calculates partial g(r) for any element pair across multi-frame trajectories
- **Angle distribution function** — sorts bond angles by element triplet and plots histograms with kernel-density-estimation fits
- **Coordination heatmap** — generates Seaborn heatmaps showing the frequency of element combinations in each coordination type
- **VASP → XYZ converter** — reads POSCAR / CONTCAR files (Direct or Cartesian) and converts them to XYZ, auto-filling cell dimensions
- **Embedded plots** — all figures render inside the GUI with full Matplotlib navigation (zoom, pan, save)
- **Chinese / English UI** — one-click language toggle in the top-right corner

## Quick Start

### Prerequisites

- Python 3.9 or later
- tkinter (included with standard Python installers; see platform notes below)

### Setup

**macOS / Linux**

```bash
git clone https://github.com/OldJii/BELLO.git
cd BELLO
chmod +x setup.sh
./setup.sh
```

**Windows**

```powershell
git clone https://github.com/OldJii/BELLO.git
cd BELLO
setup.bat
```

The setup script creates a virtual environment, installs all dependencies, and checks for tkinter availability.

### Run

```bash
# macOS / Linux
source .venv/bin/activate
python BELLO_GUI.py

# Windows
.venv\Scripts\activate.bat
python BELLO_GUI.py
```

### tkinter Troubleshooting

| Platform | Command |
|----------|---------|
| macOS (Homebrew) | `brew install python-tk@3.XX` (match your Python version) |
| Ubuntu / Debian | `sudo apt install python3-tk` |
| Fedora | `sudo dnf install python3-tkinter` |
| Windows | Re-run the Python installer → check **tcl/tk and IDLE** |

## Usage Guide

### 1. BELLO Analysis

1. Click **Browse** and select an XYZ trajectory file.
2. Enter the unit-cell dimensions (X, Y, Z in Angstroms).
3. Set the inter-atomic distance **Threshold** and **Tolerance**, or check **Automatic threshold** to let the program determine them from the RDF first peak.
4. Select **BELLO Analysis** and click **Calculate**.
5. When finished, the right panel shows:
   - **Plot 1** — Order-parameter *q* distribution (smoothed, per fold type)
   - **Plot 2** — Local-order population vs. frame number
6. Output files are written to the working directory (`output-*.txt`, `out2.pdb`, etc.).

### 2. Radial Distribution Function (RDF)

1. Load an XYZ file and fill in cell dimensions as above.
2. Select **RDF** — a dialog asks for maximum radius, step size, and the two element symbols.
3. Click **Calculate**. The partial g(r) curve appears in the plot panel; data is saved to `RDF.txt`.

### 3. Angle Distribution

> Requires a prior BELLO run (needs `output-angle-distribution.txt`).

1. Select **Angle Distribution** — enter the element symbols present in your system.
2. Click **Calculate**. Per-triplet histograms with KDE fits appear as separate tabs. Sorted data is saved under `Sorted_angles/`.

### 4. Coordination Heatmap

> Requires a prior BELLO run (needs `output-human-readable-coords.txt`).

1. Select **Coordination Heatmap** — enter element symbols.
2. Click **Calculate**. Heatmaps show the frequency of each element combination in 3-fold through octahedral environments.

### 5. VASP → XYZ Conversion

1. In the **VASP → XYZ Converter** card, click **Browse** and select a `.vasp` / `POSCAR` / `CONTCAR` file.
2. Click **Convert**. The converted `.xyz` file path is auto-filled into the input field, and cell dimensions are populated from the lattice vectors.

### 6. Language Switching

Click the **中文** / **English** button in the top-right corner. The entire interface rebuilds in the selected language instantly.

## Project Structure

```
BELLO/
├── BELLO_GUI.py                          # GUI (tkinter, cross-platform)
├── BELLO_main.py                         # Core analysis engine
├── Radial_Pair_Distribution_Function.py  # Partial RDF calculation
├── Angle_Distribution_Function.py        # Angle sorting and plotting
├── Coordination_Heatmap.py               # Coordination heatmap generation
├── vasp_converter.py                     # VASP POSCAR/CONTCAR → XYZ
├── test.xyz                              # Sample trajectory (270 Ge atoms)
├── requirements.txt                      # Python dependencies
├── setup.sh                              # macOS / Linux setup script
├── setup.bat                             # Windows setup script
└── .gitignore
```

## Dependencies

All installed automatically by the setup script:

- numpy
- pandas
- scipy
- matplotlib
- seaborn

"""
Robust XYZ trajectory file reader.

Replaces pd.read_fwf for XYZ parsing — handles multi-element files,
arbitrary comment lines, and multi-frame trajectories reliably.
"""

import numpy as np


def read_xyz(filepath):
    """Read an XYZ trajectory file and return an array compatible with BELLO internals.

    Returns a 2D object array with 4 columns and the same row layout
    that pd.read_fwf + fillna('x') would produce for a well-formed file:
      - Count lines:   [natom_str, 'x', 'x', 'x']
      - Comment lines: [comment_str, 'x', 'x', 'x']
      - Data lines:    [element_str, x_float, y_float, z_float]
    """
    with open(filepath, 'r') as fh:
        raw_lines = fh.readlines()

    rows = []
    i = 0
    while i < len(raw_lines):
        line = raw_lines[i].strip()
        if not line:
            i += 1
            continue

        try:
            natom = int(line)
        except ValueError:
            i += 1
            continue

        rows.append([str(natom), 'x', 'x', 'x'])

        i += 1
        if i < len(raw_lines):
            comment = raw_lines[i].strip() or 'x'
            rows.append([comment, 'x', 'x', 'x'])
        i += 1

        for _ in range(natom):
            if i >= len(raw_lines):
                break
            parts = raw_lines[i].split()
            i += 1
            if len(parts) >= 4:
                rows.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])
            else:
                rows.append([raw_lines[i - 1].strip(), 'x', 'x', 'x'])

    return np.array(rows, dtype=object)

"""
VASP file to XYZ converter.

Supported formats:
  - POSCAR / CONTCAR (single-frame, VASP 5+ with element symbols)
  - XDATCAR (multi-frame trajectory, VASP 5+)

Handles both 'Direct' (fractional) and 'Cartesian' coordinate types.
"""

import numpy as np
import os


def _parse_header(lines):
    """Parse the common VASP header (lines 0-6) shared by POSCAR and XDATCAR.

    Returns (comment, lattice, element_types, counts, total_atoms, header_end_idx).
    header_end_idx points to the line *after* the counts line (i.e. line 7).
    """
    comment = lines[0]
    scaling = float(lines[1])

    lattice = np.zeros((3, 3))
    for i in range(3):
        lattice[i] = [float(x) for x in lines[2 + i].split()]
    lattice *= scaling

    tokens_line5 = lines[5].split()
    try:
        int(tokens_line5[0])
        raise ValueError(
            "VASP 4 format (no element line) is not supported. "
            "Please add element symbols on line 6 of your POSCAR."
        )
    except ValueError as e:
        if "not supported" in str(e):
            raise
        element_types = tokens_line5

    counts = [int(x) for x in lines[6].split()]
    if len(element_types) != len(counts):
        raise ValueError(
            f"Element count mismatch: {len(element_types)} element names "
            f"but {len(counts)} counts."
        )

    total_atoms = sum(counts)
    element_list = []
    for elem, count in zip(element_types, counts):
        element_list.extend([elem] * count)

    return comment, lattice, element_types, counts, total_atoms, element_list, 7


def _detect_format(filepath):
    """Return 'xdatcar' or 'poscar' by inspecting file content."""
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    _, _, _, _, total_atoms, _, hdr_end = _parse_header(lines)

    idx = hdr_end
    if idx >= len(lines):
        return 'poscar'

    first_after_header = lines[idx].strip()
    if first_after_header.lower().startswith('direct configuration') or \
       first_after_header.lower().startswith('direct configuration'):
        return 'xdatcar'

    remaining_data_lines = len(lines) - hdr_end
    if first_after_header[0].upper() == 'S':
        remaining_data_lines -= 1
        idx += 1
    if idx < len(lines) and lines[idx][0].upper() in ('D', 'C', 'K'):
        remaining_data_lines -= 1

    if remaining_data_lines > total_atoms * 1.5:
        for line in lines[hdr_end:]:
            if line.lower().startswith('direct configuration') or \
               line.lower().startswith('direct config'):
                return 'xdatcar'

    return 'poscar'


def read_poscar(filepath):
    """Parse a VASP POSCAR/CONTCAR file and return structured data."""
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    comment, lattice, element_types, counts, total_atoms, element_list, hdr_end = \
        _parse_header(lines)

    next_line_idx = hdr_end
    selective_dynamics = False
    if lines[next_line_idx][0].upper() == 'S':
        selective_dynamics = True
        next_line_idx += 1

    coord_type = lines[next_line_idx][0].upper()
    if coord_type not in ('D', 'C', 'K'):
        raise ValueError(f"Unknown coordinate type: {lines[next_line_idx]}")
    is_cartesian = coord_type in ('C', 'K')

    positions = np.zeros((total_atoms, 3))
    for i in range(total_atoms):
        parts = lines[next_line_idx + 1 + i].split()
        positions[i] = [float(parts[j]) for j in range(3)]

    if not is_cartesian:
        positions = positions @ lattice

    return {
        'comment': comment,
        'lattice': lattice,
        'elements': element_list,
        'counts': counts,
        'element_types': element_types,
        'positions': positions,
        'selective_dynamics': selective_dynamics,
        'total_atoms': total_atoms,
    }


def read_xdatcar(filepath):
    """Parse a VASP XDATCAR file and return all frames.

    Returns dict with keys: comment, lattice, element_types, counts,
    elements, total_atoms, frames (list of Nx3 arrays in Cartesian).
    """
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

    comment, lattice, element_types, counts, total_atoms, element_list, hdr_end = \
        _parse_header(lines)

    frames = []
    idx = hdr_end
    while idx < len(lines):
        tag = lines[idx].lower()
        if tag.startswith('direct') or tag.startswith('d'):
            idx += 1
            if idx + total_atoms > len(lines):
                break
            positions = np.zeros((total_atoms, 3))
            for i in range(total_atoms):
                parts = lines[idx + i].split()
                positions[i] = [float(parts[j]) for j in range(3)]
            positions = positions @ lattice
            frames.append(positions)
            idx += total_atoms
        else:
            idx += 1

    if not frames:
        raise ValueError(
            "No frames found in XDATCAR. "
            "Expected 'Direct configuration= N' lines."
        )

    return {
        'comment': comment,
        'lattice': lattice,
        'element_types': element_types,
        'counts': counts,
        'elements': element_list,
        'total_atoms': total_atoms,
        'frames': frames,
    }


def _make_output_path(input_path, suffix='.xyz'):
    """Generate output path, handling files without extensions (POSCAR, CONTCAR, XDATCAR)."""
    base, ext = os.path.splitext(input_path)
    if ext:
        return base + suffix
    return input_path + suffix


def poscar_to_xyz(input_path, output_path=None):
    """Convert a VASP POSCAR/CONTCAR file to XYZ format.

    Returns:
        Tuple of (output_path, cell_x, cell_y, cell_z).
    """
    data = read_poscar(input_path)

    if output_path is None:
        output_path = _make_output_path(input_path)

    natoms = data['total_atoms']
    elements = data['elements']
    positions = data['positions']
    lattice = data['lattice']

    cell_x = np.linalg.norm(lattice[0])
    cell_y = np.linalg.norm(lattice[1])
    cell_z = np.linalg.norm(lattice[2])

    with open(output_path, 'w') as f:
        f.write(f"  {natoms}\n")
        f.write(f" {data['comment']}\n")
        for i in range(natoms):
            f.write(
                f"{elements[i]:2s}  {positions[i][0]:14.6f}"
                f"  {positions[i][1]:14.6f}"
                f"  {positions[i][2]:14.6f}\n"
            )

    return output_path, cell_x, cell_y, cell_z


def xdatcar_to_xyz(input_path, output_path=None):
    """Convert a VASP XDATCAR (multi-frame) to a multi-frame XYZ trajectory.

    Returns:
        Tuple of (output_path, cell_x, cell_y, cell_z, num_frames).
    """
    data = read_xdatcar(input_path)

    if output_path is None:
        output_path = _make_output_path(input_path)

    natoms = data['total_atoms']
    elements = data['elements']
    lattice = data['lattice']
    frames = data['frames']

    cell_x = np.linalg.norm(lattice[0])
    cell_y = np.linalg.norm(lattice[1])
    cell_z = np.linalg.norm(lattice[2])

    with open(output_path, 'w') as f:
        for frame_idx, positions in enumerate(frames):
            f.write(f"  {natoms}\n")
            f.write(f" Frame {frame_idx + 1} | {data['comment']}\n")
            for i in range(natoms):
                f.write(
                    f"{elements[i]:2s}  {positions[i][0]:14.6f}"
                    f"  {positions[i][1]:14.6f}"
                    f"  {positions[i][2]:14.6f}\n"
                )

    return output_path, cell_x, cell_y, cell_z, len(frames)


def convert_vasp(input_path, output_path=None):
    """Auto-detect VASP file type and convert to XYZ.

    Returns:
        Tuple of (output_path, cell_x, cell_y, cell_z, num_frames).
        num_frames is 1 for POSCAR/CONTCAR.
    """
    fmt = _detect_format(input_path)
    if fmt == 'xdatcar':
        return xdatcar_to_xyz(input_path, output_path)
    else:
        path, cx, cy, cz = poscar_to_xyz(input_path, output_path)
        return path, cx, cy, cz, 1

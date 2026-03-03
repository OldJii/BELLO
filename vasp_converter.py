"""
VASP POSCAR/CONTCAR (.vasp) to XYZ file converter.

Supports VASP 5+ format with element symbols on line 6.
Handles both 'Direct' (fractional) and 'Cartesian' coordinate types.
"""

import numpy as np
import os


def read_poscar(filepath):
    """Parse a VASP POSCAR/CONTCAR file and return structured data.

    Returns dict with keys: comment, scaling, lattice, elements,
    counts, coord_type, positions, selective_dynamics.
    """
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines()]

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
        elements = tokens_line5

    counts = [int(x) for x in lines[6].split()]
    if len(elements) != len(counts):
        raise ValueError(
            f"Element count mismatch: {len(elements)} element names "
            f"but {len(counts)} counts."
        )

    next_line_idx = 7
    selective_dynamics = False
    if lines[next_line_idx][0].upper() == 'S':
        selective_dynamics = True
        next_line_idx = 8

    coord_type = lines[next_line_idx][0].upper()
    if coord_type not in ('D', 'C', 'K'):
        raise ValueError(f"Unknown coordinate type: {lines[next_line_idx]}")
    is_cartesian = coord_type in ('C', 'K')

    total_atoms = sum(counts)
    positions = np.zeros((total_atoms, 3))
    for i in range(total_atoms):
        parts = lines[next_line_idx + 1 + i].split()
        positions[i] = [float(parts[j]) for j in range(3)]

    if not is_cartesian:
        positions = positions @ lattice

    element_list = []
    for elem, count in zip(elements, counts):
        element_list.extend([elem] * count)

    return {
        'comment': comment,
        'lattice': lattice,
        'elements': element_list,
        'counts': counts,
        'element_types': elements,
        'positions': positions,
        'selective_dynamics': selective_dynamics,
        'total_atoms': total_atoms,
    }


def poscar_to_xyz(input_path, output_path=None):
    """Convert a VASP POSCAR/CONTCAR file to XYZ format.

    Args:
        input_path: Path to the .vasp / POSCAR / CONTCAR file.
        output_path: Path for the output .xyz file.
                     If None, replaces extension with .xyz.

    Returns:
        Tuple of (output_path, cell_x, cell_y, cell_z) where cell dimensions
        are the diagonal elements of the lattice matrix (orthorhombic approx).
    """
    data = read_poscar(input_path)

    if output_path is None:
        base, _ = os.path.splitext(input_path)
        output_path = base + '.xyz'

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


def batch_poscar_to_xyz(input_paths, output_dir=None):
    """Convert multiple VASP files to XYZ, returning list of results."""
    results = []
    for path in input_paths:
        out = None
        if output_dir:
            base = os.path.splitext(os.path.basename(path))[0]
            out = os.path.join(output_dir, base + '.xyz')
        results.append(poscar_to_xyz(path, out))
    return results

import shutil
import numpy as np
import os
import pathlib
import stat
import glob
import math as mt
import matplotlib.pyplot as plt
import scipy.stats as stats
from xyz_reader import read_xyz


def folder_overwrite(sorted_folder):
    if os.path.exists(sorted_folder):
        def _remove_readonly(func, path, _exc_info):
            os.chmod(path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
            if func is os.open:
                func(path, os.O_RDONLY)
            else:
                func(path)
        shutil.rmtree(sorted_folder, onerror=_remove_readonly)
    os.makedirs(sorted_folder, exist_ok=True)


def compute_adf(url, celldmx, celldmy, celldmz,
                center_elem, ligand_elem_1, ligand_elem_2,
                cutoff, progress_callback=None,
                frame_stride=1, max_frames=0):
    """Compute angle distribution for a specific triplet directly from XYZ.

    Finds all angles ligand1-center-ligand2 where both bonds are within *cutoff*.

    Args:
        url: path to XYZ file
        celldmx/y/z: cell dimensions
        center_elem: center atom element symbol (e.g. 'Ge')
        ligand_elem_1: first ligand element (e.g. 'Se')
        ligand_elem_2: second ligand element (e.g. 'Se')
        cutoff: bond distance cutoff in Å
        progress_callback: optional callable(stage, cur, total)
        frame_stride, max_frames: trajectory sampling

    Returns:
        Tuple of (angles_array, fig) where angles_array is 1-D numpy array
        of all bond angles in degrees, and fig is a matplotlib Figure.
    """
    def _pcb(stage, cur, total):
        if progress_callback:
            progress_callback(stage, cur, total)

    frame_stride = int(frame_stride)
    max_frames = int(max_frames)
    max_frames_arg = None if max_frames == 0 else max_frames
    file = read_xyz(url, frame_stride=frame_stride, max_frames=max_frames_arg)
    if len(file) == 0:
        raise ValueError("No valid frames were parsed from XYZ file.")
    Natom = int(file[0][0])
    lfile = len(file)
    Nframes = lfile // (Natom + 2)
    if Nframes <= 0:
        raise ValueError("No frames available after parsing/sampling.")

    cx, cy, cz = float(celldmx), float(celldmy), float(celldmz)
    cutoff = float(cutoff)

    fa_rows = []
    for i in range(lfile):
        if file[i, 1] != 'x':
            x, y, z = file[i, 1], file[i, 2], file[i, 3]
            if x > cx or x < 0: x = x % cx
            if y > cy or y < 0: y = y % cy
            if z > cz or z < 0: z = z % cz
            fa_rows.append([file[i, 0], x, y, z])
    fa = np.array(fa_rows, dtype=object)
    if len(fa) == 0:
        raise ValueError("No atom coordinates found in XYZ data.")

    def pbc_distance(a, b):
        dx = abs(a[0] - b[0]); dx = min(dx, cx - dx)
        dy = abs(a[1] - b[1]); dy = min(dy, cy - dy)
        dz = abs(a[2] - b[2]); dz = min(dz, cz - dz)
        return mt.sqrt(dx**2 + dy**2 + dz**2)

    def pbc_vector(center, neighbor):
        """Return the minimum-image vector from center to neighbor."""
        v = np.array([float(neighbor[j] - center[j]) for j in range(3)])
        cell = [cx, cy, cz]
        for j in range(3):
            if v[j] > cell[j] / 2:
                v[j] -= cell[j]
            elif v[j] < -cell[j] / 2:
                v[j] += cell[j]
        return v

    all_angles = []
    triplet_name = f"{ligand_elem_1}-{center_elem}-{ligand_elem_2}"
    sym_same = (ligand_elem_1 == ligand_elem_2)

    for N in range(0, len(fa), Natom):
        _pcb(f'ADF {triplet_name}', N, len(fa))
        fo = fa[N:N + Natom]
        lfo = len(fo)

        for c_idx in range(lfo):
            if fo[c_idx, 0] != center_elem:
                continue
            center_pos = fo[c_idx, 1:4].astype(float)

            neighbors_1 = []
            neighbors_2 = []
            for n_idx in range(lfo):
                if n_idx == c_idx:
                    continue
                d = pbc_distance(center_pos, fo[n_idx, 1:4].astype(float))
                if d > cutoff or d < 1e-8:
                    continue
                if fo[n_idx, 0] == ligand_elem_1:
                    neighbors_1.append((n_idx, d))
                if fo[n_idx, 0] == ligand_elem_2:
                    neighbors_2.append((n_idx, d))

            for i, (n1_idx, d1) in enumerate(neighbors_1):
                start_j = (i + 1) if sym_same else 0
                for j in range(start_j, len(neighbors_2)):
                    n2_idx, d2 = neighbors_2[j]
                    if n1_idx == n2_idx:
                        continue
                    v1 = pbc_vector(center_pos, fo[n1_idx, 1:4].astype(float))
                    v2 = pbc_vector(center_pos, fo[n2_idx, 1:4].astype(float))
                    cos_angle = np.dot(v1, v2) / (d1 * d2)
                    cos_angle = max(-1.0, min(1.0, cos_angle))
                    angle_deg = mt.degrees(mt.acos(cos_angle))
                    all_angles.append(angle_deg)

    angles = np.array(all_angles)
    fig = _plot_adf(angles, triplet_name)
    return angles, fig


def _plot_adf(angles, triplet_name):
    """Create a histogram + KDE plot for an angle distribution."""
    fig = plt.figure(facecolor="none", figsize=(6.75, 5))
    ax = fig.add_subplot(111)
    fontsize = 14
    ax.set_ylabel('Distribution (arb. units)', size=fontsize)
    ax.set_xlabel('Angle (°)', size=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize)

    if angles.size == 0:
        ax.text(0.5, 0.5, f'No angles found for {triplet_name}',
                transform=ax.transAxes, ha='center', va='center', fontsize=12)
        fig.tight_layout()
        return fig

    bins = np.linspace(0, 180, 45, dtype='i')
    n, x, patches = ax.hist(x=angles, bins=bins, label=triplet_name,
                            edgecolor='black', density=True)
    cm = plt.colormaps['viridis']
    n_range = n.max() - n.min()
    if n_range > 0:
        col = (n - n.min()) / n_range
    else:
        col = np.zeros_like(n)
    for c_val, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm(c_val))

    if angles.size > 1:
        density = stats.gaussian_kde(angles)
        ax.plot(x, density(x), label=f'KDE {triplet_name}')
    else:
        ax.axvline(float(angles[0]), linestyle='--', linewidth=2,
                   label=f'Single sample {triplet_name}')

    ax.legend()
    fig.tight_layout()
    return fig


def sorter(elements):
    """Sort angle distributions by element triplets from BELLO output and plot.

    Returns:
        List of matplotlib Figure objects.
    """
    plt.rcParams['figure.max_open_warning'] = 0

    file = np.loadtxt('output-angle-distribution.txt', delimiter=' ',
                       dtype=str, usecols=(0, 1), ndmin=2)
    length = file.shape[0]

    sorted_folder = pathlib.Path(os.getcwd(), 'Sorted_angles')
    folder_overwrite(sorted_folder)

    names = []

    def func_element(elements):
        sorted_list = []
        trash_collection = []
        condition = False
        for i in elements:
            for j in elements:
                for k in elements:
                    name1 = i + '-' + j + '-' + k
                    name2 = k + '-' + j + '-' + i
                    condition = False
                    for x in range(0, length):
                        if (file[x][1] == name1 or file[x][1] == name2) and (file[x][1] not in trash_collection):
                            sorted_list.append(list(file[x]))
                            temp_list = [i, k]
                            temp_list.sort()
                            name = temp_list[0] + '-' + j + '-' + temp_list[1]
                            condition = True
                    trash_collection.append(name1)
                    trash_collection.append(name2)
                    if condition:
                        save_name = name + '.txt'
                        names.append(name)
                        np.savetxt(str(sorted_folder / save_name),
                                   sorted_list, fmt="%s")
                        sorted_list = []
                        condition = False

    func_element(elements)
    print('Angle sorting is done!')
    bins = np.linspace(0, 180, 45, dtype='i')

    file_names = glob.glob(str(sorted_folder / '*.txt'))
    file_names = [os.path.basename(fn) for fn in file_names]

    figures = []
    for i in range(0, len(file_names)):
        fig = plt.figure(facecolor="none", figsize=(6.75, 5))
        ax = fig.add_subplot(111)
        ax.set_ylabel('Distribution (arb. units)', size=14)
        ax.set_xlabel('Angles (Theta)', size=14)

        data = np.loadtxt(str(sorted_folder / file_names[i]),
                          dtype="f", usecols=0)
        data = np.atleast_1d(data).astype(float)
        if data.size == 0:
            plt.close(fig)
            continue
        n, x, patches = ax.hist(x=data, bins=bins, label=names[i],
                                edgecolor='black', density=True)
        cm = plt.colormaps['viridis']
        n_range = n.max() - n.min()
        if n_range > 0:
            col = (n - n.min()) / n_range
        else:
            col = np.zeros_like(n)
        for c_val, p in zip(col, patches):
            plt.setp(p, 'facecolor', cm(c_val))
        fontsize = 14
        ax.tick_params(axis='both', labelsize=fontsize)
        if data.size > 1:
            density = stats.gaussian_kde(data)
            ax.plot(x, density(x), label=('Fit %s' % names[i]))
        else:
            ax.axvline(float(data[0]), linestyle='--', linewidth=2,
                       label=('Single sample %s' % names[i]))
        plt.legend()
        plt.tight_layout()
        figures.append(fig)

    return figures

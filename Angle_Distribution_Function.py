import shutil
import numpy as np
import os
import pathlib
import stat
import glob
import matplotlib.pyplot as plt
import scipy.stats as stats


def folder_overwrite(sorted_folder):
    if os.path.exists(sorted_folder):
        def _remove_readonly(func, path, _exc_info):
            os.chmod(path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
            func(path)
        shutil.rmtree(sorted_folder, onerror=_remove_readonly)
    os.makedirs(sorted_folder, exist_ok=True)


def sorter(elements):
    """Sort angle distributions by element triplets and plot histograms.

    Returns:
        List of matplotlib Figure objects.
    """
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

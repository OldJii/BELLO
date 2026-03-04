import numpy as np
import pandas as pd
from collections import Counter
from itertools import combinations_with_replacement as cr
import seaborn as sns
import matplotlib.pyplot as plt


def coordination_heatmap(*elements_raw):
	"""Build coordination heatmaps from BELLO output.

	Returns:
		List of matplotlib Figure objects.
	"""
	rows = []
	with open('output-human-readable-coords.txt', 'r') as fh:
		for line in fh:
			parts = line.split()
			if len(parts) >= 2:
				rows.append((parts[0], parts[1]))
	if not rows:
		raise ValueError("output-human-readable-coords.txt is empty or invalid.")
	file = np.array(rows, dtype='str')
	Nframe = int(file[0, 1])

	def merger(input_list):
		if isinstance(input_list, list):
			input_list.sort()
		else:
			input_list = sorted(input_list)
		output = ''.join([f'{key}{value}' for key, value in Counter(input_list).items()])
		return output

	def string_attach(input_list):
		output = []
		for i in input_list:
			output.append(merger(i))
		return output

	def column_length(elements, fold):
		return len(list(cr(elements, fold)))

	def indexer(grid):
		fold3_index = string_attach(list(cr(elements, 3)))
		grid.loc[:column_length(elements, 3), '3_fold_index'] = fold3_index
		fold4_index = string_attach(list(cr(elements, 4)))
		grid.loc[:column_length(elements, 4), '4_fold_index'] = fold4_index
		grid.loc[:column_length(elements, 4), 'Tetrahedral_index'] = fold4_index
		fold5_index = string_attach(list(cr(elements, 5)))
		grid.loc[:column_length(elements, 5), '5_fold_index'] = fold5_index
		fold6_index = string_attach(list(cr(elements, 6)))
		grid.loc[:column_length(elements, 6), '6_fold_index'] = fold6_index
		return grid

	def heatmap_params(num_rows):
		"""Choose a readable figure size/font based on combination count."""
		height = max(6.0, min(26.0, 2.5 + 0.28 * float(num_rows)))
		if num_rows > 60:
			annot_size = 5
		elif num_rows > 35:
			annot_size = 6
		elif num_rows > 20:
			annot_size = 7
		else:
			annot_size = 8
		return (10.0, height), annot_size

	def draw_heatmap(data_grid, title):
		fig_size, annot_size = heatmap_params(len(data_grid.index))
		fig, ax = plt.subplots(figsize=fig_size)
		sns.heatmap(
			data_grid.loc[:, ('3_fold', '4_fold', 'Tetrahedral', '5_fold', '6_fold')],
			annot=data_grid.loc[:, ('3_fold_index', '4_fold_index', 'Tetrahedral_index', '5_fold_index', '6_fold_index')],
			linewidths=.5,
			cmap="viridis",
			fmt='',
			square=False,
			yticklabels=False,
			annot_kws={'size': annot_size},
			ax=ax,
		)
		ax.set_title(title)
		fig.tight_layout()
		return fig

	elements = []
	for item in elements_raw:
		if isinstance(item, (list, tuple)):
			elements.extend(item)
		else:
			elements.append(item)
	elements.sort()

	columns = ['3_fold_index', '3_fold', '4_fold_index', '4_fold',
	           'Tetrahedral_index', 'Tetrahedral', '5_fold_index', '5_fold',
	           '6_fold_index', '6_fold']
	index_cols = [c for c in columns if c.endswith('_index')]
	count_cols = [c for c in columns if not c.endswith('_index')]
	length = column_length(elements, 6)
	idx = np.linspace(1, length, length, dtype='i')
	grid = pd.DataFrame(0.0, index=idx, columns=count_cols)
	for col in index_cols:
		grid[col] = ''
	grid = grid.reindex(columns=columns)
	grid = indexer(grid)

	figures = []

	for i in range(0, len(file)):
		if file[i, 0] == '3-FOLD':
			temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0]]
			for j in range(1, column_length(elements, 3)):
				if merger(temp_element) == grid.loc[j, '3_fold_index']:
					grid.loc[j, '3_fold'] += 1
		elif file[i, 0] == '4-FOLD':
			temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0]]
			for j in range(1, column_length(elements, 4)):
				if merger(temp_element) == grid.loc[j, '4_fold_index']:
					grid.loc[j, '4_fold'] += 1
		elif file[i, 0] == 'TETRAHEDRAL':
			temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0]]
			for j in range(1, column_length(elements, 4)):
				if merger(temp_element) == grid.loc[j, 'Tetrahedral_index']:
					grid.loc[j, 'Tetrahedral'] += 1
		elif file[i, 0] == '5-FOLD':
			temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0], file[i + 6, 0]]
			for j in range(1, column_length(elements, 5)):
				if merger(temp_element) == grid.loc[j, '5_fold_index']:
					grid.loc[j, '5_fold'] += 1
		elif file[i, 0] == 'OCTAHEDRAL':
			temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0], file[i + 6, 0], file[i + 7, 0]]
			for j in range(1, column_length(elements, 6)):
				if merger(temp_element) == grid.loc[j, '6_fold_index']:
					grid.loc[j, '6_fold'] += 1

	grid.loc[:, ('3_fold', '4_fold', 'Tetrahedral', '5_fold', '6_fold')] = (
		grid.loc[:, ('3_fold', '4_fold', 'Tetrahedral', '5_fold', '6_fold')]) / Nframe

	fig_global = draw_heatmap(grid, "Coordination Heatmap (all center atoms)")
	figures.append(fig_global)

	for x in elements:
		grid = pd.DataFrame(0.0, index=idx, columns=count_cols)
		for col in index_cols:
			grid[col] = ''
		grid = grid.reindex(columns=columns)
		grid = indexer(grid)
		for i in range(0, len(file)):
			if file[i, 0] == '3-FOLD' and file[i + 1, 0] == x:
				temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0]]
				for j in range(1, column_length(elements, 3)):
					if merger(temp_element) == grid.loc[j, '3_fold_index']:
						grid.loc[j, '3_fold'] += 1
			elif file[i, 0] == '4-FOLD' and file[i + 1, 0] == x:
				temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0]]
				for j in range(1, column_length(elements, 4)):
					if merger(temp_element) == grid.loc[j, '4_fold_index']:
						grid.loc[j, '4_fold'] += 1
			elif file[i, 0] == 'TETRAHEDRAL' and file[i + 1, 0] == x:
				temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0]]
				for j in range(1, column_length(elements, 4)):
					if merger(temp_element) == grid.loc[j, 'Tetrahedral_index']:
						grid.loc[j, 'Tetrahedral'] += 1
			elif file[i, 0] == '5-FOLD' and file[i + 1, 0] == x:
				temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0], file[i + 6, 0]]
				for j in range(1, column_length(elements, 5)):
					if merger(temp_element) == grid.loc[j, '5_fold_index']:
						grid.loc[j, '5_fold'] += 1
			elif file[i, 0] == 'OCTAHEDRAL' and file[i + 1, 0] == x:
				temp_element = [file[i + 2, 0], file[i + 3, 0], file[i + 4, 0], file[i + 5, 0], file[i + 6, 0], file[i + 7, 0]]
				for j in range(1, column_length(elements, 6)):
					if merger(temp_element) == grid.loc[j, '6_fold_index']:
						grid.loc[j, '6_fold'] += 1
		grid.loc[:, ('3_fold', '4_fold', 'Tetrahedral', '5_fold', '6_fold')] = (
			grid.loc[:, ('3_fold', '4_fold', 'Tetrahedral', '5_fold', '6_fold')]) / Nframe
		fig_elem = draw_heatmap(grid, "Coordination Heatmap with %s as center atom" % x)
		figures.append(fig_elem)

	return figures

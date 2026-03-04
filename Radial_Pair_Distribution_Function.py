import numpy as np
import math as mt
import matplotlib.pyplot as plt
from xyz_reader import read_xyz


def RDF(url, celldmx_raw, celldmy_raw, celldmz_raw, dr, rmax,
        first_element, second_element, progress_callback=None,
        frame_stride=1, max_frames=0):
	"""Compute partial radial pair distribution function g(r).

	Returns:
		List containing one matplotlib Figure.
	"""
	def _progress(stage, cur, total):
		if progress_callback:
			progress_callback(stage, cur, total)

	print("|-----------------------------------------------------|\n"
	      "|----------------------B.E.L.L.O----------------------|\n"
	      "|---------Bond Element Lattice Locality Order---------|\n"
	      "|-----------------------------------------------------|\n"
	      "|----------Radial Pair Distribution Function----------|")

	frame_stride = int(frame_stride)
	max_frames = int(max_frames)
	if frame_stride < 1:
		raise ValueError(f"frame_stride must be >= 1, got {frame_stride}")
	if max_frames < 0:
		raise ValueError(f"max_frames must be >= 0, got {max_frames}")

	max_frames_arg = None if max_frames == 0 else max_frames
	file = read_xyz(url, frame_stride=frame_stride, max_frames=max_frames_arg)
	if len(file) == 0:
		raise ValueError("No valid frames were parsed from XYZ file.")
	Natom = int(file[0][0])

	celldmx = float(celldmx_raw)
	celldmy = float(celldmy_raw)
	celldmz = float(celldmz_raw)
	dr = float(dr)
	rmax = float(rmax)

	fat = first_element
	sat = second_element

	lfile = len(file)
	print("File length is: ", lfile)

	Nframes = lfile // (Natom + 2)
	print("Number of frames: ", Nframes)
	if Nframes <= 0:
		raise ValueError("No frames available after parsing/sampling.")

	fa_rows = []

	for i in range(lfile):
		if file[i, 1] != 'x':
			x = file[i, 1]
			y = file[i, 2]
			z = file[i, 3]
			if x > celldmx or x < 0:
				x = x % celldmx
			if y > celldmy or y < 0:
				y = y % celldmy
			if z > celldmz or z < 0:
				z = z % celldmz
			fa_rows.append([file[i, 0], x, y, z])
	fa = np.array(fa_rows, dtype=object)
	print("Boundary Condition is done!\nCalculating Radial Pair Distribution Function:")
	if len(fa) == 0:
		raise ValueError("No atom coordinates found in XYZ data.")

	def distance(a, b):
		dx = abs(a[0] - b[0])
		x = min(dx, abs(A - dx))
		dy = abs(a[1] - b[1])
		y = min(dy, abs(B - dy))
		dz = abs(a[2] - b[2])
		z = min(dz, abs(C - dz))
		return mt.sqrt(x**2 + y**2 + z**2)

	A = celldmx
	B = celldmy
	C = celldmz
	totalRDF = []
	n = int(mt.ceil(rmax / dr))
	r = np.around(np.arange(0 + dr, rmax + dr, dr), decimals=3)
	vol = [(4.0 / 3.0) * mt.pi * (r[i])**3 for i in range(n)]
	vol.insert(0, 0)

	for N in range(0, len(fa), Natom):
		_progress('Frames', N, len(fa))
		rlist = []
		natm1 = 0
		natm2 = 0
		fo = np.copy(fa[N:N + Natom])
		lfo = len(fo)

		for l in range(lfo):
			if fo[l, 0] == fat:
				natm1 += 1
			if fo[l, 0] == sat:
				natm2 += 1

		for c in range(0, lfo):
			if fo[c, 0] == fat:
				for b1 in range(0, lfo):
					if (c != b1) & (fo[b1, 0] == sat):
						a = fo[c, 1:4]
						b = fo[b1, 1:4]
						db1 = distance(a, b)
						if (db1 <= rmax) and (db1 <= celldmx / 2) and (db1 <= celldmy / 2) and (db1 <= celldmz / 2):
							rlist.append(db1)

		if not rlist:
			continue

		rdf_hist = np.zeros(n)
		density = natm2 / (celldmx * celldmy * celldmz)
		if density == 0 or natm1 == 0:
			continue
		hist_counts, _ = np.histogram(rlist, bins=np.arange(0, rmax + dr, dr))
		rdf_hist[:len(hist_counts)] = hist_counts.astype(float)

		for i in range(n):
			rdf_hist[i] /= natm1
			rdf_hist[i] /= density
			rdf_hist[i] /= vol[i] - vol[i - 1]
		totalRDF.append(rdf_hist)

	if not totalRDF:
		raise ValueError(f"No RDF pairs found for element pair {fat}-{sat}.")

	rdf = np.mean(np.array(totalRDF), axis=0)
	finalfile = np.column_stack((r, rdf))
	np.savetxt('RDF.txt', finalfile, delimiter=' ', fmt="%s")

	fig = plt.figure(facecolor="none", figsize=(6.75, 5))
	ax = fig.add_subplot(111)
	fontsize = 14
	ax.tick_params(axis='both', labelsize=fontsize)
	plt.plot(r, rdf, label='Radial distribution function', lw=2)
	plt.xlabel('Radius (Angstroms)', fontsize=fontsize)
	plt.ylabel('Radial pair distribution function g(r)', fontsize=fontsize)
	plt.legend()
	plt.tight_layout()

	return [fig]

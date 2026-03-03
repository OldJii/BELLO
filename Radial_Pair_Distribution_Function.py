import numpy as np
import math as mt
import pandas as pd
import matplotlib.pyplot as plt


def RDF(url, celldmx_raw, celldmy_raw, celldmz_raw, dr, rmax,
        first_element, second_element, progress_callback=None):
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

	f = pd.read_fwf(url, header=None)
	f = f.fillna("x")
	file = np.array(f)
	Natom = int(f[0][0])

	celldmx = float(celldmx_raw)
	celldmy = float(celldmy_raw)
	celldmz = float(celldmz_raw)
	dr = float(dr)
	rmax = float(rmax)

	fat = first_element
	sat = second_element

	lfile = len(file)
	print("File length is: ", lfile)

	fa = np.empty((0, 4), int)
	templist = []
	intervals = 1 / dr
	Nframes = round(len(f) / (Natom + 2))
	print("Number of frames: ", Nframes)

	for i in range(lfile):
		if file[i, 1] != 'x':
			if file[i, 1] > celldmx or file[i, 1] < 0:
				file[i, 1] = file[i, 1] % celldmx
			if file[i, 2] > celldmy or file[i, 2] < 0:
				file[i, 2] = file[i, 2] % celldmy
			if file[i, 3] > celldmz or file[i, 3] < 0:
				file[i, 3] = file[i, 3] % celldmz
			fa = np.vstack((fa, file[i]))
	print("Boundary Condition is done!\nCalculating Radial Pair Distribution Function:")

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
						templist.append(db1)
						if (db1 <= rmax) and (db1 <= celldmx / 2) and (db1 <= celldmy / 2) and (db1 <= celldmz / 2):
							rlist.append(db1)

		if not rlist:
			continue

		rmax_frame = max(rlist)
		nvals = len(rlist)
		n = int(mt.ceil(rmax_frame / dr))
		rdf_hist = np.zeros(n)
		rr = np.arange(0 + dr, rmax_frame + dr, dr)
		r = np.around(rr, decimals=3)
		vol = [(4.0 / 3.0) * mt.pi * (r[i])**3 for i in range(n)]
		vol.insert(0, 0)
		density = natm2 / (celldmx * celldmy * celldmz)

		for i in range(n):
			for j in range(nvals):
				test = round(rlist[j] * intervals) / intervals
				if test == r[i]:
					rdf_hist[i] += 1

		for i in range(n):
			rdf_hist[i] /= natm1
			rdf_hist[i] /= density
			rdf_hist[i] /= vol[i] - vol[i - 1]
		totalRDF.append(rdf_hist)

	rdf = sum(totalRDF) / Nframes
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

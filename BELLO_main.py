import numpy as np
import math as mt
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import Akima1DInterpolator


def BELLO(url, celldmx_raw, celldmy_raw, celldmz_raw, automatic,
          trhld_raw, tlrnc_raw, angle_condition, progress_callback=None):
	"""Run BELLO local-order analysis.

	Args:
		progress_callback: Optional callable(stage_name, current, total)
			for reporting progress to the GUI.

	Returns:
		List of matplotlib Figure objects produced by the analysis.
	"""
	def _progress(stage, cur, total):
		if progress_callback:
			progress_callback(stage, cur, total)

	print("|-----------------------------------------------------|\n"
	      "|----------------------B.E.L.L.O----------------------|\n"
	      "|---------Bond Element Lattice Locality Order---------|\n"
	      "|-----------------------------------------------------|")

	f = pd.read_fwf(url, header=None)
	f = f.fillna("x")
	file = np.array(f)
	Natom = int(f[0][0])
	lfile = len(file)
	celldmx = float(celldmx_raw)
	celldmy = float(celldmy_raw)
	celldmz = float(celldmz_raw)

	print("File length is: ", lfile)

	fa = np.empty((0, 4), int)
	templist = []
	Nframes = round(len(f) / (Natom + 2))
	print("Number of frames: ", round(Nframes))

	def gaus(x, a, x0, sigma):
		return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

	for i in range(lfile):
		_progress('Boundary condition', i, lfile)
		if file[i, 1] != 'x':
			if file[i, 1] > celldmx or file[i, 1] < 0:
				file[i, 1] = file[i, 1] % celldmx
			if file[i, 2] > celldmy or file[i, 2] < 0:
				file[i, 2] = file[i, 2] % celldmy
			if file[i, 3] > celldmz or file[i, 3] < 0:
				file[i, 3] = file[i, 3] % celldmz
			fa = np.vstack((fa, file[i]))
	print("Boundary Condition is done!")

	if round(len(fa) - Natom) == 0:
		max_val = 1
	else:
		max_val = round(len(fa))

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

	figures = []

	def automatic_rdf(fa, Natom):
		for N in range(0, len(fa), Natom):
			_progress('Automatic Threshold', N, max_val)
			rlist = []
			fo = np.copy(fa[N:N + Natom])
			lfo = len(fo)
			dr = 0.1
			rmax = 10
			intervals = 1 / dr
			natm1 = Natom
			natm2 = Natom

			for c in range(0, lfo):
				for b1 in range(0, lfo):
					a = fo[c, 1:4]
					b = fo[b1, 1:4]
					db1 = distance(a, b)
					templist.append(db1)
					if (db1 <= rmax) and (db1 <= celldmx / 2) and (db1 <= celldmy / 2) and (db1 <= celldmz / 2):
						rlist.append(db1)

			rmax = max(rlist)
			nvals = len(rlist)
			n = int(mt.ceil(rmax / dr))
			rdf_hist = np.zeros(n)
			rr = np.arange(0 + dr, rmax + dr, dr)
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

			max_rdf = rdf.argmax()
			g_range1 = max(0, max_rdf - 5)
			g_range2 = min(len(r), max_rdf + 5)

			gx = r[g_range1:g_range2]
			gy = rdf[g_range1:g_range2]

			x = len(gx)
			mean = 3
			sigma = sum(gy * (gx - mean)**2) / x

			popt, pcov = curve_fit(gaus, gx, gy, p0=[1, mean, sigma])
			print('Peak value: ', popt[0])
			print('Gaussian peak location: ', popt[1])
			print('Sigma: ', popt[2])

			fig_rdf = plt.figure(facecolor="none", figsize=(6.75, 5))
			plt.plot(r, rdf_hist, label='Radial distribution function')
			plt.plot(gx, gaus(gx, *popt), 'ro:', label='Gaussian fit')
			plt.xlabel('Radius (Angstroms)')
			plt.ylabel('Radial pair distribution function g(r)')
			plt.legend()
			plt.tight_layout()
			figures.append(fig_rdf)

			trhld = popt[1]
			if popt[2] * 0.7 < 0.1:
				tlrnc = 0.1
			else:
				tlrnc = popt[2] * 0.7
			return trhld, tlrnc

	if automatic:
		trhld, tlrnc = automatic_rdf(fa, Natom)
	else:
		trhld = float(trhld_raw)
		tlrnc = float(tlrnc_raw)

	atomnumbers = []
	coordinates = []
	coordinates.append(['Total_frames', Nframes, ' ', ' '])
	mcoordinates = []
	pdbcoords = []
	xyzcoordinates = []
	localstats = []
	localstats.append(["0-FOLD", "1-FOLD", "2-FOLD", "3-FOLD", "4-FOLD",
	                    "TETRAHEDRAL", "5-FOLD", "OCTAHEDRAL", "Total"])
	t1PDB = "PDB coordinates"
	t2PDB = "-------BELLO-------"
	t3PDB = "CRYST1   %f   %f   %f  90.00  90.00  90.00 P1          1" % (celldmx, celldmy, celldmz)
	pdbcoords.append(t1PDB)
	pdbcoords.append(t2PDB)
	pdbcoords.append(t3PDB)
	qframe = []
	q3flframe = []
	q4flframe = []
	qtetframe = []
	q5flframe = []
	q6flframe = []
	tetatot = []

	def rep(a):
		binn = []
		binn.append(a[0])
		for i in range(1, len(a)):
			if np.any([a[i][1:4] != x[1:4] for x in binn]):
				binn.append(a[i])
		return binn

	def dott(a, b, c):
		celld = [celldmx, celldmy, celldmz]
		tlist = [a, b]
		for j in tlist:
			d = mt.sqrt((j[0] - c[0])**2 + (j[1] - c[1])**2 + (j[2] - c[2])**2)
			if d != distance(j, c):
				for i in range(0, 3):
					if abs(j[i] - c[i]) > (celld[i] / 2):
						if c[i] < j[i]:
							j[i] = j[i] - celld[i]
						elif c[i] > j[i]:
							j[i] = j[i] + celld[i]
		return sum((tlist[0] - c) * (tlist[1] - c))

	for N in range(0, len(fa), Natom):
		_progress('Frames', N, max_val)
		fo = np.copy(fa[N:N + Natom])
		lfo = len(fo)
		q = []
		q3fl = []
		q4fl = []
		qtet = []
		q5fl = []
		q6fl = []
		aq = []
		mcoordinates = []
		nlocal = 0
		n0fl = 0
		n1fl = 0
		n2fl = 0
		n3fl = 0
		n4fl = 0
		ntet = 0
		n5fl = 0
		noct = 0
		for c in range(0, lfo):
			templist = []
			condition = 'false'
			nbond = 0
			if condition == 'false':
				for b1 in range(0, lfo):
					if c != b1 and condition == 'false':
						a = np.copy(fo[c, 1:4])
						b = np.copy(fo[b1, 1:4])
						ab1 = np.copy(fo[b1, 1:4])
						db1 = distance(a, b)
						if (trhld - tlrnc) < db1 < (trhld + tlrnc):
							nbond += 1
							if condition == 'false':
								for b2 in range(0, lfo):
									if c != b2 and b1 != b2 and condition == 'false':
										b = np.copy(fo[b2, 1:4])
										ab2 = np.copy(fo[b2, 1:4])
										db2 = distance(a, b)
										if abs(db1 - db2) < tlrnc:
											nbond += 1
											dot = dott(ab1, ab2, a)
											cosb12 = dot / (db1 * db2)
											cosb12 = max(-1.0, min(1.0, cosb12))
											teta = mt.degrees(mt.acos(cosb12))

											if condition == 'false':
												for b3 in range(0, lfo):
													if b3 != c and b3 != b2 and b3 != b1 and condition == 'false':
														b = np.copy(fo[b3, 1:4])
														ab3 = np.copy(fo[b3, 1:4])
														db3 = distance(a, b)
														if abs(db3 - db1) <= tlrnc:
															nbond += 1
															tetatot.append(
																"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b1, 0],
																                                   fo[c, 0], fo[b2, 0]))
															dot = dott(ab3, ab2, a)
															cosb23 = dot / (db3 * db2)
															cosb23 = max(-1.0, min(1.0, cosb23))
															teta = mt.degrees(mt.acos(cosb23))
															tetatot.append(
																"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b2, 0],
																                                   fo[c, 0], fo[b3, 0]))
															dot = dott(ab1, ab3, a)
															cosb13 = dot / (db3 * db1)
															cosb13 = max(-1.0, min(1.0, cosb13))
															teta = mt.degrees(mt.acos(cosb13))
															tetatot.append(
																"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b1, 0],
																                                   fo[c, 0], fo[b3, 0]))

															if condition == 'false':
																for b4 in range(0, lfo):
																	if c != b4 and b1 != b4 and b2 != b4 and b3 != b4 and condition == 'false':
																		b = np.copy(fo[b4, 1:4])
																		ab4 = np.copy(fo[b4, 1:4])
																		db4 = distance(a, b)
																		if abs(db4 - db1) <= tlrnc:
																			nbond += 1
																			dot = dott(ab1, ab4, a)
																			cosb14 = dot / (db1 * db4)
																			cosb14 = max(-1.0, min(1.0, cosb14))
																			teta = mt.degrees(mt.acos(cosb14))
																			tetatot.append(
																				"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b1, 0], fo[c, 0], fo[b4, 0]))
																			dot = dott(ab2, ab4, a)
																			cosb24 = dot / (db2 * db4)
																			cosb24 = max(-1.0, min(1.0, cosb24))
																			teta = mt.degrees(mt.acos(cosb24))
																			tetatot.append(
																				"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b2, 0], fo[c, 0], fo[b4, 0]))
																			dot = dott(ab4, ab3, a)
																			cosb34 = dot / (db3 * db4)
																			cosb34 = max(-1.0, min(1.0, cosb34))
																			teta = mt.degrees(mt.acos(cosb34))
																			tetatot.append(
																				"{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b3, 0], fo[c, 0], fo[b4, 0]))

																			if condition == 'false':
																				for b5 in range(0, lfo):
																					if c != b5 and b1 != b5 and b2 != b5 and b3 != b5 and b4 != b5 and condition == 'false':
																						b = np.copy(fo[b5, 1:4])
																						ab5 = np.copy(fo[b5, 1:4])
																						db5 = distance(a, b)
																						if abs(db5 - db1) <= tlrnc:
																							nbond += 1
																							dot = dott(ab1, ab5, a)
																							cosb15 = dot / (db1 * db5)
																							cosb15 = max(-1.0, min(1.0, cosb15))
																							teta = mt.degrees(mt.acos(cosb15))
																							tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b1, 0], fo[c, 0], fo[b5, 0]))
																							dot = dott(ab2, ab5, a)
																							cosb25 = dot / (db2 * db5)
																							cosb25 = max(-1.0, min(1.0, cosb25))
																							teta = mt.degrees(mt.acos(cosb25))
																							tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b2, 0], fo[c, 0], fo[b5, 0]))
																							dot = dott(ab3, ab5, a)
																							cosb35 = dot / (db3 * db5)
																							cosb35 = max(-1.0, min(1.0, cosb35))
																							teta = mt.degrees(mt.acos(cosb35))
																							tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b3, 0], fo[c, 0], fo[b5, 0]))
																							dot = dott(ab4, ab5, a)
																							cosb45 = dot / (db4 * db5)
																							cosb45 = max(-1.0, min(1.0, cosb45))
																							teta = mt.degrees(mt.acos(cosb45))
																							tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b4, 0], fo[c, 0], fo[b5, 0]))

																							if condition == 'false':
																								for b6 in range(0, lfo):
																									if c != b6 and b1 != b6 and b2 != b6 and b3 != b6 and b4 != b6 and b5 != b6 and condition == 'false':
																										b = np.copy(fo[b6, 1:4])
																										ab6 = np.copy(fo[b6, 1:4])
																										db6 = distance(a, b)
																										if abs(db6 - db1) <= tlrnc:
																											nbond += 1
																											dot = dott(ab1, ab6, a)
																											cosb16 = dot / (db1 * db6)
																											cosb16 = max(-1.0, min(1.0, cosb16))
																											teta = mt.degrees(mt.acos(cosb16))
																											tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b1, 0], fo[c, 0], fo[b6, 0]))
																											dot = dott(ab2, ab6, a)
																											cosb26 = dot / (db2 * db6)
																											cosb26 = max(-1.0, min(1.0, cosb26))
																											teta = mt.degrees(mt.acos(cosb26))
																											tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b2, 0], fo[c, 0], fo[b6, 0]))
																											dot = dott(ab3, ab6, a)
																											cosb36 = dot / (db3 * db6)
																											cosb36 = max(-1.0, min(1.0, cosb36))
																											teta = mt.degrees(mt.acos(cosb36))
																											tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b3, 0], fo[c, 0], fo[b6, 0]))
																											dot = dott(ab4, ab6, a)
																											cosb46 = dot / (db4 * db6)
																											cosb46 = max(-1.0, min(1.0, cosb46))
																											teta = mt.degrees(mt.acos(cosb46))
																											tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b4, 0], fo[c, 0], fo[b6, 0]))
																											dot = dott(ab5, ab6, a)
																											cosb56 = dot / (db5 * db6)
																											cosb56 = max(-1.0, min(1.0, cosb56))
																											teta = mt.degrees(mt.acos(cosb56))
																											tetatot.append("{:3.3f} {:2s}-{:2s}-{:2s}".format(teta, fo[b5, 0], fo[c, 0], fo[b6, 0]))
																											condition = 'true'
																										else:
																											continue
																									if condition == 'true':
																										break
																								condition = 'true'
																								break
																							else:
																								continue
																						else:
																							continue
																					if condition == 'true':
																						break
																				condition = 'true'
																				break
																			else:
																				continue
																		else:
																			continue
																	if condition == 'true':
																		break
																condition = 'true'
																break
														else:
															continue
													if condition == 'true':
														break
												condition = 'true'
												break
											else:
												continue
										else:
											continue
									if condition == 'true':
										break
								condition = 'true'
								break
							else:
								continue
						else:
							continue
					if condition == 'true':
						break
				else:
					condition = 'true'
				if nbond == 0:
					nlocal += 1
					n0fl += 1
					q = "Not defined"
					templist = [c]
					coordinates.append(["0-FOLD", " -", "q is:", q])
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('0-FOLD', c + 1, 0, 0, 0, 0, 0, 0)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], '0FL', 'A', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 1:
					nlocal += 1
					n1fl += 1
					q = "Not defined"
					templist = [c, b1]
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('1-FOLD', c + 1, b1 + 1, 0, 0, 0, 0, 0)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], '1FL', 'B', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 2:
					nlocal += 1
					n2fl += 1
					q = "Not defined"
					templist = [c, b1, b2]
					coordinates.append(["2-FOLD", " -", "q is:", q])
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('2-FOLD', c + 1, b1 + 1, b2 + 1, 0, 0, 0, 0)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], '2FL', 'C', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 3:
					nlocal += 1
					n3fl += 1
					q = 1 - (0.375 * ((cosb12 + 1 / 3)**2 + (cosb13 + 1 / 3)**2 + (cosb23 + 1 / 3)**2))
					q3fl.append(q)
					aq.append(q)
					templist = [c, b1, b2, b3]
					coordinates.append(["3-FOLD", " -", "q is:", q])
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('3-FOLD', c + 1, b1 + 1, b2 + 1, b3 + 1, 0, 0, 0)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], '3FL', 'D', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 4:
					q = 1 - (0.375 * ((cosb12 + 1 / 3)**2 + (cosb13 + 1 / 3)**2 + (cosb14 + 1 / 3)**2 +
					                   (cosb23 + 1 / 3)**2 + (cosb24 + 1 / 3)**2 + (cosb34 + 1 / 3)**2))
					aq.append(q)
					templist = [c, b1, b2, b3, b4]
					if 1 >= q >= 0.85:
						qtet.append(q)
						nlocal += 1
						ntet += 1
						coordinates.append(["TETRAHEDRAL", " -", "q is:", q])
						atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('TETHDL', c + 1, b1 + 1, b2 + 1, b3 + 1, b4 + 1, 0, 0)
						atomnumbers.append(atomnumbertemp)
						for i in templist:
							mcoordinates.append(fo[i, 0:4])
							coordinates.append(fo[i, 0:4])
							pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
								'ATOM', i + 1, fo[i, 0], 'TET', 'E', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
							pdbcoords.append(pdbtemp)
						coordinates.append([" ", " ", " ", " "])
					elif q <= 0.85:
						q4fl.append(q)
						nlocal += 1
						n4fl += 1
						coordinates.append(["4-FOLD", " -", "q is:", q])
						atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('4-FOLD', c + 1, b1 + 1, b2 + 1, b3 + 1, b4 + 1, 0, 0)
						atomnumbers.append(atomnumbertemp)
						for i in templist:
							mcoordinates.append(fo[i, 0:4])
							coordinates.append(fo[i, 0:4])
							pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
								'ATOM', i + 1, fo[i, 0], '4FL', 'F', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
							pdbcoords.append(pdbtemp)
						coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 5:
					nlocal += 1
					n5fl += 1
					q = 1 - (0.375 * ((cosb12 + 1 / 3)**2 + (cosb13 + 1 / 3)**2 + (cosb14 + 1 / 3)**2 +
					                   (cosb23 + 1 / 3)**2 + (cosb24 + 1 / 3)**2 + (cosb34 + 1 / 3)**2 +
					                   (cosb15 + 1 / 3)**2 + (cosb25 + 1 / 3)**2 + (cosb35 + 1 / 3)**2 +
					                   (cosb45 + 1 / 3)**2))
					q5fl.append(q)
					aq.append(q)
					templist = [c, b1, b2, b3, b4, b5]
					coordinates.append(["5-FOLD", " -", "q is:", q])
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('5-FOLD', c + 1, b1 + 1, b2 + 1, b3 + 1, b4 + 1, b5 + 1, 0)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], '5FL', 'G', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
				if nbond == 6:
					nlocal += 1
					noct += 1
					q = 1 - (0.375 * ((cosb12 + 1 / 3)**2 + (cosb13 + 1 / 3)**2 + (cosb14 + 1 / 3)**2 +
					                   (cosb23 + 1 / 3)**2 + (cosb24 + 1 / 3)**2 + (cosb34 + 1 / 3)**2 +
					                   (cosb15 + 1 / 3)**2 + (cosb25 + 1 / 3)**2 + (cosb35 + 1 / 3)**2 +
					                   (cosb45 + 1 / 3)**2 + (cosb16 + 1 / 3)**2 + (cosb26 + 1 / 3)**2 +
					                   (cosb36 + 1 / 3)**2 + (cosb46 + 1 / 3)**2 + (cosb56 + 1 / 3)**2))
					q6fl.append(q)
					aq.append(q)
					templist = [c, b1, b2, b3, b4, b5, b6]
					coordinates.append(["OCTAHEDRAL", " -", "q is:", q])
					atomnumbertemp = '{:6s} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d} {:03d}'.format('OCTHDL', c + 1, b1 + 1, b2 + 1, b3 + 1, b4 + 1, b5 + 1, b6 + 1)
					atomnumbers.append(atomnumbertemp)
					for i in templist:
						mcoordinates.append(fo[i, 0:4])
						coordinates.append(fo[i, 0:4])
						pdbtemp = '{:6s}{:5d} {:^4s} {:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(
							'ATOM', i + 1, fo[i, 0], 'OCT', 'H', nlocal, fo[i, 1], fo[i, 2], fo[i, 3], 1.00, 0.00, fo[i, 0], "00")
						pdbcoords.append(pdbtemp)
					coordinates.append([" ", " ", " ", " "])
					condition = 'true'
			else:
				continue
		Endf = "End of frame %d " % (N / Natom)
		if len(mcoordinates) != 0:
			for x in rep(mcoordinates):
				xyzcoordinates.append(x)
		xyzcoordinates.append(["\r\n", Endf, "\r\n", " "])
		atomnumbers.append(Endf)
		stat = " \r\n Number of 0-FOLD, 1-FOLD, 2-FOLD, 3-FOLD, 4-FOLD, TETRAHEDRAL, 5-FOLD, OCTAHEDRAL: %d, %d, %d, %d, %d, %d, %d, %d" % (
			n0fl, n1fl, n2fl, n3fl, n4fl, ntet, n5fl, noct)
		locq = "\r\n Number of local orders are: %d " % nlocal
		coordinates.append([Endf, stat, locq, "\r\n"])
		localstats.append([n0fl, n1fl, n2fl, n3fl, n4fl, ntet, n5fl, noct, nlocal])
		pdbcoords.append(Endf)
		aq = [round(num, 2) for num in aq]
		qframe.append(aq)
		q3fl = [round(num, 2) for num in q3fl]
		q3flframe.append(q3fl)
		q4fl = [round(num, 2) for num in q4fl]
		q4flframe.append(q4fl)
		qtet = [round(num, 2) for num in qtet]
		qtetframe.append(qtet)
		q5fl = [round(num, 2) for num in q5fl]
		q5flframe.append(q5fl)
		q6fl = [round(num, 2) for num in q6fl]
		q6flframe.append(q6fl)
		print(Endf)

	if angle_condition:
		cond = (N / Natom) % 500
		if (cond == 0) or ((N / Natom) == int((len(fa) / Natom) - 1)):
			txtt = "Angle_distribution_frame_%d.txt" % (N / Natom)
			np.savetxt(txtt, tetatot, fmt="%s")
			tetatot = []

	np.savetxt('output-q3fl.txt', np.array(q3flframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-q4fl.txt', np.array(q4flframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-qtet.txt', np.array(qtetframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-q5fl.txt', np.array(q5flframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-q6fl.txt', np.array(q6flframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-q-total.txt', np.array(qframe, dtype=object), delimiter=',', fmt="%s")
	np.savetxt('output-human-readable-coords.txt', coordinates, delimiter=' ', fmt="%s")
	np.savetxt('output-atom-number.txt', atomnumbers, delimiter=' ', fmt='%s')
	np.savetxt('output-xyz-coords.txt', mcoordinates, delimiter=' ', fmt='%s')
	np.savetxt('output-pdb-coords-m.txt', pdbcoords, fmt="%s")
	np.savetxt('out2.pdb', pdbcoords, fmt="%s")
	np.savetxt('output-local-statistics.txt', localstats, fmt="%s")
	np.savetxt('output-angle-distribution.txt', tetatot, fmt="%s")
	print("Done!")

	bin_number = 20
	initial_range = 0

	q_histogram_values_tot = np.concatenate(qframe, axis=None)
	q_histogram_tot, bins_tot = np.histogram(a=q_histogram_values_tot, range=(initial_range, 1.0), bins=bin_number)
	bins_tot = list(np.linspace(initial_range, 1, len(bins_tot) - 1))
	x_data_smooth = np.linspace(initial_range, 1, 1000)
	y_non_smooth_tot = Akima1DInterpolator(bins_tot, q_histogram_tot)
	y_smooth_tot = y_non_smooth_tot(x_data_smooth)

	q_histogram_values_3fl = np.concatenate(q3flframe, axis=None)
	q_histogram_3fl, _ = np.histogram(a=q_histogram_values_3fl, range=(initial_range, 1.0), bins=bin_number)
	y_smooth_3fl = Akima1DInterpolator(bins_tot, q_histogram_3fl)(x_data_smooth)

	q_histogram_values_4fl = np.concatenate(q4flframe, axis=None)
	q_histogram_4fl, _ = np.histogram(a=q_histogram_values_4fl, range=(initial_range, 1.0), bins=bin_number)
	y_smooth_4fl = Akima1DInterpolator(bins_tot, q_histogram_4fl)(x_data_smooth)

	q_histogram_values_tet = np.concatenate(qtetframe, axis=None)
	q_histogram_tet, _ = np.histogram(a=q_histogram_values_tet, range=(initial_range, 1.0), bins=bin_number)
	y_smooth_tet = Akima1DInterpolator(bins_tot, q_histogram_tet)(x_data_smooth)

	q_histogram_values_5fl = np.concatenate(q5flframe, axis=None)
	q_histogram_5fl, _ = np.histogram(a=q_histogram_values_5fl, range=(initial_range, 1.0), bins=bin_number)
	y_smooth_5fl = Akima1DInterpolator(bins_tot, q_histogram_5fl)(x_data_smooth)

	q_histogram_values_6fl = np.concatenate(q6flframe, axis=None)
	q_histogram_6fl, _ = np.histogram(a=q_histogram_values_6fl, range=(initial_range, 1.0), bins=bin_number)
	y_smooth_6fl = Akima1DInterpolator(bins_tot, q_histogram_6fl)(x_data_smooth)

	edgecolor = (0, 0, 0, 1)

	fig_q = plt.figure(facecolor="none", figsize=(6.75, 5))
	ax = fig_q.add_subplot(111)
	text_label_position = (5 / 6) * max(q_histogram_tot) if max(q_histogram_tot) > 0 else 1
	fontsize = 14
	fontsize_legend = 10
	ax.tick_params(axis='x', labelsize=fontsize)
	ax.tick_params(axis='y', labelleft=False)

	plt.fill_between(x_data_smooth, y_smooth_tot, edgecolor=edgecolor, facecolor=(1, 1, 1, 1), label='Total', lw=2)
	plt.fill_between(x_data_smooth, y_smooth_3fl, edgecolor=edgecolor, alpha=0.5, label='3_Fold', lw=2)
	plt.fill_between(x_data_smooth, y_smooth_4fl, edgecolor=edgecolor, alpha=0.5, label='4_Fold', lw=2)
	plt.fill_between(x_data_smooth, y_smooth_tet, edgecolor=edgecolor, alpha=0.5, label='Tetrahedral', lw=2)
	plt.fill_between(x_data_smooth, y_smooth_5fl, edgecolor=edgecolor, alpha=0.5, label='5_Fold', lw=2)
	plt.fill_between(x_data_smooth, y_smooth_6fl, edgecolor=edgecolor, alpha=0.5, label='Octahedral', lw=2)

	ref_lines = [
		(0.75, '3FL planar'), (0.875, '3FL pyramidal'), (0.625, '4FL defective'),
		(0.5, '4FL planar'), (1.0, 'Tetrahedral'), (0.333, '5FL defective'), (0.0, 'Octahedral'),
	]
	for xval, label in ref_lines:
		plt.axvline(xval, color='grey', ls='--', lw=2)
		plt.text(xval, text_label_position, label, color='grey', rotation=90,
		         size=fontsize_legend, ha='center', va='center',
		         bbox=dict(facecolor='w', edgecolor='w', boxstyle='round'))

	plt.ylabel("Distribution (arb. units)", size=fontsize)
	plt.xlabel('Order parameter q', size=fontsize)
	plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
	           mode="expand", borderaxespad=0, ncol=3)
	plt.tight_layout()
	figures.append(fig_q)

	labels = ["0-FOLD", "1-FOLD", "2-FOLD", "3-FOLD", "4-FOLD", "TETRAHDRL", "5-FOLD", "6-FOLD"]
	y_axis = np.array(localstats[1:])
	y_axis = y_axis[:, :-1].astype(float)
	x_axis = list(np.linspace(0, Nframes, Nframes))

	fig_stats = plt.figure(facecolor="none", figsize=(6.75, 5))
	for i in range(len(labels)):
		plt.plot(x_axis, y_axis[:, i], label=labels[i], marker='o')
	plt.yticks(list(np.linspace(0, np.max(y_axis), 10)), fontsize=14)
	plt.xticks(fontsize=14)
	plt.xlabel('Frame number', size=fontsize)
	plt.ylabel('Local order population / frame', size=fontsize)
	plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
	           mode="expand", borderaxespad=0, ncol=4)
	plt.grid(visible=True, axis='y')
	plt.tight_layout()
	figures.append(fig_stats)

	return figures


import numpy as np 
import matplotlib.pyplot as plt
import random

def dist(p1, p2):
	d = 0
	for i in range(3):
		d += (p1[i] - p2[i]) ** 2
	d = np.sqrt(d)
	return d

def check_for_overlap(p1, p2, r):
	if dist(p1, p2) < 2 * r:
		return True
	else:
		return False

def convert_to_cart(r, theta, phi):
	x = r * np.cos(theta) * np.sin(phi)
	y = r * np.sin(theta) * np.sin(phi)
	z = r * np.cos(phi)
	return x,y,z

def make_grid_coords(nmol, d):
	coords = []
	if nmol == 16:
		rect = [[d/2, d/2], [-d/2, d/2], [d/2, -d/2], [-d/2, -d/2]]
		for z in [-d, 0, d, 2*d]:
			for pt in rect:
				coord = (pt[0], pt[1], z)
				coords.append(coord)

	elif nmol == 27:
		rect = [
		[0, 0], [d, 0], [0, d], 
		[d, d], [-d, 0], [-d, -d], 
		[0, -d], [-d, d], [d, -d]
		]
		for z in [-d, 0, d]:
			for pt in rect:
				coord = (pt[0], pt[1], z)
				coords.append(coord)

	elif nmol == 125:
		x = -3*d
		y = -3*d
		rect = []
		for i in range(5):
			x += d
			y = -3*d
			for j in range(5):
				y += d
				rect.append([x, y])
		for z in [-2*d, -d, 0, d, 2*d]:
			for pt in rect:
				coord = (pt[0], pt[1], z)
				coords.append(coord)
	return coords

def make_rand_coords(nmol, molr, maxr):
	coords = []
	i = 0
	while i < nmol:
		r = maxr*np.random.rand()
		theta = 2*np.pi*np.random.rand()
		phi = np.pi*np.random.rand()
		test = convert_to_cart(r, theta, phi)
		if any([check_for_overlap(coord, test, molr) for coord in coords]):
			continue
		else:
			coords.append(test)
			i += 1
	return coords

def write_to_xyz(nmol, coords, outname=None):
	if outname is None:
		outname = "pos%i_grid.xyz" % nmol
	f = open(outname, 'w+')
	for coord in coords:
		x, y, z = coord
		line = "%f\t%f\t%f\n" % (x, y, z)
		f.write(line)
	f.close()

def write_3bd_infile(nmol, salt):
	descriptor = "_%i_%.2f" % (nmol, salt)
	f = open("run.manybodyapprox%s.inp" % descriptor, 'w+')
	f.write("runname bodyapprox%s\n" % descriptor )
	f.write("runtype bodyapprox\n")
	f.write("3bdloc 3bd%s.dat\n" % descriptor)
	f.write("2bdloc 2bd%s.dat\n" % descriptor)
	f.write("salt %.2f\n" % salt)
	f.write("temp 298\n")
	f.write("idiel 4\n")
	f.write("sdiel 78\n")
	# f.write("randorient\n")
	f.write("attypes 1\n")
	f.write("units kT\n")
	f.write("type 1 %i\n" % nmol)
	f.write("pqr 1 mol.pqr\n")
	f.write("xyz 1 pos%i_grid.xyz\n" % nmol)
	f.close()

def write_energyforce_infile(nmol, salt, mix=False):
	descriptor = "_%i_%.2f" % (nmol, salt)
	f = open("run.energyforce%s.inp" % descriptor, 'w+')
	f.write("runtype energyforce 1\n")
	f.write("runname energyforce%s.out\n" % descriptor )
	f.write("units kT\n")
	f.write("salt %.2f\n" % salt)
	f.write("temp 298\n")
	f.write("idiel 4\n")
	f.write("sdiel 78\n")
	# f.write("randorient\n")
	f.write("attypes 1\n")
	f.write("type 1 %i\n" % nmol)
	f.write("pqr 1 mol_pos.pqr\n")
	f.write("xyz 1 pos%i_grid.xyz\n" % nmol)
	f.close()

def split_coords(coords, m=2, rand=False):
	csets = [[] for _ in range(m)]
	if rand:
		random.shuffle(coords)

	i = 0
	while i < len(coords):
		for j in range(m):
			csets[j].append(coords[i])
			i += 1
			if i >= len(coords):
				break
	return csets


num_mols = [16, 27, 125]
salts = [0.00, 0.05, 0.10]

mol_rad = 2
max_rad = 20

for n in num_mols:
	# rand_coords = make_rand_coords(n, mol_rad, max_rad)
	grid_coords = make_grid_coords(n, 5)
	split = split_coords(grid_coords)
	write_to_xyz(n, split[0], outname="pos%i_grid_half1.xyz" % n)
	write_to_xyz(n, split[1], outname="pos%i_grid_half2.xyz" % n)
	# for s in salts:
	# 	write_3bd_infile(n, s)
	# 	write_energyforce_infile(n, s)

# energy_force_timings = [(1.360217, 1.360518, 1.374350), 
# 						(8.289805, 8.606490, 9.198424), 
# 						(2048.688232, 1997.770020, 2024.487305)]


# threebod_timings = [(8.066987, 8.221177, 8.060463),
# 					(41.614868, 41.433807, 41.466129),
# 					(4118.210938, 4168.622559, 4173.109863)]


# def make_time_plot(nmols, enforce_times, threebod_times):
# 	fig = plt.figure()
# 	ax = fig.add_subplot(111)
# 	x = []
# 	y1 = []
# 	y2 = []
# 	for i in range(len(nmols)):
# 		av_en_time = sum(enforce_times[i]) / 3
# 		av_3bd_time = sum(threebod_times[i]) / 3
# 		x.append(nmols[i])
# 		y1.append(av_en_time)
# 		y2.append(av_3bd_time)

# 	ax.plot(x, np.log(y1), color='b', label="Complete")
# 	ax.plot(x, np.log(y2), color='g', label="Three Body approximation")
# 	ax.scatter(x, np.log(y1), s=70, color='b')
# 	ax.scatter(x, np.log(y2), s=70, color='g')

# 	ax.set_xticks(num_mols)
# 	ax.set_xlabel("Number of molecules", fontsize=16)
# 	ax.set_ylabel("Log of Calculation time (s)", fontsize=16)
# 	ax.legend(frameon=False, fontsize=16, loc='upper left')
# 	plt.tight_layout()
# 	plt.savefig("timings.pdf")

# make_time_plot(num_mols, energy_force_timings, threebod_timings)




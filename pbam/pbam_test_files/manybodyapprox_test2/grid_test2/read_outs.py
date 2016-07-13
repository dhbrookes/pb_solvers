import matplotlib.pyplot as plt
import numpy as np

def read_en_force_out(fname):
	f = open(fname)
	lines = f.readlines()
	f.close()
	ens = []
	fs = []
	tors = []
	for line in lines:
		if 'ENERGY' in line:
			split = line.split()
			ens.append(float(split[1]))
		elif 'FORCE' in line:
			split = line.split()
			fs.append(float(split[1].strip(',')))
		elif 'TORQUE' in line:
			split = line.split()
			tors.append(float(split[1].strip(',')))
	return ens, fs, tors

def read_3bd_out(fname):
	f = open(fname)
	lines = f.readlines()
	f.close()
	ens = []
	fs = []
	for line in lines:
		if 'Energy' in line:
			split = line.split()
			ens.append(float(split[1]))
			fs.append(float(split[3]))
	return ens, fs

num_mols = [16, 27, 125]
salts = [0.00, 0.05, 0.10]

# e_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]
# f_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]

en_errors = []
f_errors = []

for k in range(len(num_mols)):
	en_errors.append([])
	f_errors.append([])
	nmol = num_mols[k]
	folder = "%i_grid" % nmol
	for j in range(len(salts)):
		en_errors[k].append([])
		f_errors[k].append([])
		# if k == 2:
		# 	break
		salt = salts[j]
		descriptor = "_%i_%.2f" % (nmol, salt)
		enf_fname = "%s/energyforce%s_posneg.out" % (folder,descriptor)
		three_fname = "%s/bodyapprox%s_posneg.out" % (folder,descriptor)
		ens1, fs1, tors1 = read_en_force_out(enf_fname)
		ens2, fs2, tors2 = read_en_force_out(three_fname)

		# en_errors = []
		# f_errors = []
		for i in range(nmol):
			en_errors[k][j].append(np.abs( np.abs(ens1[i] - ens2[i]) / ens1[i] ))
			f_errors[k][j].append(np.abs( np.abs(fs1[i] - fs2[i]) / fs1[i] ))

		# e_ave_errors[j][k] = sum(en_errors) / nmol
		# f_ave_errors[j][k] = sum(f_errors) / nmol

		# print nmol, salt, sum(en_errors) / nmol, sum(f_errors) / nmol
print "Type\tN\tsalt\tmax error\taverage error\tindex of max error mols"
for i in range(len(num_mols)):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	labels = ['0.00M', '0.05M', '0.10M']
	for k in range(len(salts)):
		en_mx = max(en_errors[i][k])
		en_mx_idx = [m for m, n in enumerate(en_errors[i][k]) if n == en_mx]
		print "Energy", num_mols[i], salts[k], en_mx, sum(en_errors[i][k])/num_mols[i], [idx+1 for idx in en_mx_idx]
		# ax.hist(en_errors[i][k], bins=30, alpha=0.5, label=labels[k])
	for k in range(len(salts)):
		f_mx = max(f_errors[i][k])
		f_mx_idx = [m for m, n in enumerate(f_errors[i][k]) if n == f_mx]
		print "Force", num_mols[i], salts[k], f_mx, sum(en_errors[i][k])/num_mols[i], [idx+1 for idx in f_mx_idx]
		ax.hist(f_errors[i][k], bins=30, alpha=0.5, label=labels[k])
	print "\n"
	ax.set_ylabel('Histogram of errors')
	ax.set_xlabel('Error magnitude')
	ax.legend(loc='upper left')
	plt.savefig('hist_%i_f_posneg.png' % num_mols[i], dpi=300)
	plt.close()



# n_groups = len(num_mols)
# opacity = 1
# index = np.arange(n_groups)
# bar_width = 0.2

# rects1 = ax.bar(index, e_ave_errors[0], bar_width,
#                  alpha=opacity, label='0.00 M Salt')

# rects2 = ax.bar(index + bar_width, e_ave_errors[1], bar_width,
#                  alpha=opacity, label='0.05 M Salt', color='r')

# rects3 = ax.bar(index + bar_width*2, e_ave_errors[2], bar_width,
#                  alpha=opacity, label='0.10 M Salt', color='g')


# ax.legend(frameon=False, fontsize=14, loc='upper left')
# ax.set_xlabel("Number of molecules", fontsize=16)
# ax.set_ylabel("Average fraction error of \nthree body energy approximation", fontsize=16)
# print [idx + 1.5*bar_width for idx in index]
# ax.set_xticks([idx + 1.5*bar_width for idx in index])
# ax.set_xticklabels(num_mols)
# # ax.set_ylim(0, 0.035)
# ax.set_xlim(-0.2, 3)
# plt.tight_layout()
# plt.savefig("threebody_en_error.png", dpi=300)

# fig = plt.figure()
# ax = fig.add_subplot(111)

# n_groups = len(num_mols)
# opacity = 1
# index = np.arange(n_groups)
# bar_width = 0.2

# rects1 = ax.bar(index, f_ave_errors[0], bar_width,
#                  alpha=opacity, label='0.00 M Salt')

# rects2 = ax.bar(index + bar_width, f_ave_errors[1], bar_width,
#                  alpha=opacity, label='0.05 M Salt', color='r')

# rects3 = ax.bar(index + bar_width*2, f_ave_errors[2], bar_width,
#                  alpha=opacity, label='0.10 M Salt', color='g')


# ax.legend(frameon=False, fontsize=14, loc='upper left')
# ax.set_xlabel("Number of molecules", fontsize=16)
# ax.set_ylabel("Average fraction error of \nthree body force approximation", fontsize=16)
# ax.set_xticks([idx + 1.5*bar_width for idx in index])
# ax.set_xticklabels(num_mols)
# # ax.set_ylim(0, 0.035)
# ax.set_xlim(-0.2, 3)
# plt.tight_layout()
# plt.savefig("threebody_f_error.png", dpi=300)
# plt.close()




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

e_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]
f_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]

for k in range(len(num_mols)):
	nmol = num_mols[k]
	for j in range(len(salts)):
		if k == 2:
			break
		salt = salts[j]
		descriptor = "_%i_%.2f" % (nmol, salt)
		enf_fname = "energyforce%s.out" % descriptor 
		three_fname = "bodyapprox%s.out" % descriptor
		ens1, fs1, tors1 = read_en_force_out(enf_fname)
		ens2, fs2 = read_3bd_out(three_fname)
		print len(ens1), len(ens2), nmol
		assert len(ens1) == len(ens2) == nmol
		assert len(fs1) == len(fs2) == nmol

		en_errors = []
		f_errors = []
		for i in range(nmol):
			en_errors.append( np.abs(ens1[i] - ens2[i]) / ens1[i] )
			f_errors.append( np.abs(fs1[i] - fs2[i]) / fs1[i] )

		e_ave_errors[j][k] = sum(en_errors) / nmol
		f_ave_errors[j][k] = sum(f_errors) / nmol



		# print nmol, salt, sum(en_errors) / nmol, sum(f_errors) / nmol

fig = plt.figure()
ax = fig.add_subplot(111)

n_groups = len(num_mols)
opacity = 1
index = np.arange(n_groups)
bar_width = 0.2

rects1 = ax.bar(index, e_ave_errors[0], bar_width,
                 alpha=opacity, label='0.00 M Salt')

rects2 = ax.bar(index + bar_width, e_ave_errors[1], bar_width,
                 alpha=opacity, label='0.05 M Salt', color='r')

rects3 = ax.bar(index + bar_width*2, e_ave_errors[2], bar_width,
                 alpha=opacity, label='0.10 M Salt', color='g')


ax.legend(frameon=False, fontsize=14, loc='upper left')
ax.set_xlabel("Number of molecules", fontsize=16)
ax.set_ylabel("Average fraction error of \nthree body energy approximation", fontsize=16)
print [idx + 1.5*bar_width for idx in index]
ax.set_xticks([idx + 1.5*bar_width for idx in index])
ax.set_xticklabels(num_mols)
# ax.set_ylim(0, 0.035)
ax.set_xlim(-0.2, 2.3+0.2)
plt.tight_layout()
plt.savefig("threebody_en_error.pdf")

fig = plt.figure()
ax = fig.add_subplot(111)

n_groups = len(num_mols)
opacity = 1
index = np.arange(n_groups)
bar_width = 0.2

rects1 = ax.bar(index, f_ave_errors[0], bar_width,
                 alpha=opacity, label='0.00 M Salt')

rects2 = ax.bar(index + bar_width, f_ave_errors[1], bar_width,
                 alpha=opacity, label='0.05 M Salt', color='r')

rects3 = ax.bar(index + bar_width*2, f_ave_errors[2], bar_width,
                 alpha=opacity, label='0.10 M Salt', color='g')


ax.legend(frameon=False, fontsize=14, loc='upper left')
ax.set_xlabel("Number of molecules", fontsize=16)
ax.set_ylabel("Average fraction error of \nthree body force approximation", fontsize=16)
ax.set_xticks([idx + 1.5*bar_width for idx in index])
ax.set_xticklabels(num_mols)
# ax.set_ylim(0, 0.035)
ax.set_xlim(-0.2, 3)
plt.tight_layout()
plt.savefig("threebody_f_error.pdf")


grid_times = [(3.124841, 2.991938, 3.146087),
				(16.248241, 17.435644, 18.134485), ()
				]

three_body_times = [(13.450019, 13.489721, 13.196977), (64.520111, 66.041336, 65.543159), (5900.064453,)]




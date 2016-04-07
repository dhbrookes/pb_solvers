'''
plot trends in Rg
'''

import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/lfelberg/Desktop/test/manybodyapprox_test/")
colors = [ (1, 0.2, 0.2), (1, 0.6, 0.2), 
            (0.2, 0.6, 0.), (0.2, 0.2, 1.)]

fname = '2bd53.dat' 
ylab = fname[0]+'body energy (kT)'
dat = np.zeros((5,100000))

fig = plt.figure(figsize=(5, 3.5))
ax = fig.add_subplot(111)

ax.set_xlabel('Distance ($\AA$)', fontsize = 15)
ax.set_ylabel(ylab, fontsize = 15)
ct = 0
f = open(fname)
for line in f:
        temp = line.split()
        for i in range(len(temp)):
            if ',' in temp[i]:
                dat[i][ct] = float(temp[i][:-1])
            else:
                dat[i][ct] = float(temp[i])
        ct += 1
f.close()    

#ax.axis([2, 10, 70, 125]) 
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.xaxis.labelpad = -1
ax.yaxis.labelpad = -1

dist = dat[0][0:ct]
energy = dat[1][0:ct]
force = np.power(dat[2][0:ct]*dat[2][0:ct] 
                    + dat[3][0:ct]*dat[3][0:ct]
                    + dat[4][0:ct]*dat[4][0:ct], 0.5)
                    
if 'energy' in ylab:
    y = energy
    til = '_energy'
else:
    y = force
    til = '_force'

ax.scatter(dist, y, marker = '.',)

plt.savefig(fname[:-4] + til +'_v_dist.jpg', 
                    bbox_inches = 'tight')
plt.close()
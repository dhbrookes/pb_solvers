'''
plot trends in Rg
'''

import numpy as np
import matplotlib.pyplot as plt
import os

def getlin(fname):
  '''Function to get mbdy data'''
  ct = 0
  dat = np.zeros((5,1000000))
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
  return dat, ct


#os.chdir("/Users/felb315/Desktop/grid_test/")
os.chdir("/Users/davidbrookes/Projects/pb_solvers/pbam/pbam_test_files/manybodyapprox_test2/grid_test2/125_grid")
colors = [ 'red', 'blue', 'green', 'yellow']

mrkrs = [ 'o', 'v', 's'] 


nbdy = [2, 3]
num_mols = [125] #, 125]
salts = [0.00, 0.05, 0.10]

e_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]
f_ave_errors = [[0 for _ in range(len(num_mols))] for _ in range(len(salts))]

for i in range(len(nbdy)):
  bdy = nbdy[i]
  fig = plt.figure(figsize=(5, 3.5))
  ax = fig.add_subplot(111)
  for k in range(len(num_mols)):
    nmol = num_mols[k]
    for j in range(len(salts)):
      salt = salts[j]
      descriptor = "{0:d}bd_{1:d}_{2:.2f}.dat".format(bdy, nmol, salt)
      dat, ct = getlin(descriptor)
      dist = dat[0][0:ct]
      energy = dat[1][0:ct]
      force = np.power(dat[2][0:ct]*dat[2][0:ct] 
                    + dat[3][0:ct]*dat[3][0:ct]
                    + dat[4][0:ct]*dat[4][0:ct], 0.5)
                    
      ax.scatter(dist, force, marker = mrkrs[k], color = colors[j],
                 label="{0:d} mols, {1:.2f}M".format(nmol, salt))

  
  ylab = str(bdy)+'body energy (kT)'
  #ax.axis([2, 10, 70, 125]) 
  ax.set_xlabel('Distance ($\AA$)', fontsize = 15)
  ax.set_ylabel(ylab, fontsize = 15)
  
  ax.tick_params(axis='x', labelsize=18)
  ax.tick_params(axis='y', labelsize=18)
  ax.xaxis.labelpad = -1
  ax.yaxis.labelpad = -1
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
             ncol=2, mode="expand", borderaxespad=0.,
             fontsize=8)
  
  plt.savefig('dist_v_force' + str(bdy) +'.png', 
               dpi=300, bbox_inches = 'tight')
  plt.close()
  

# 2016 May 2 lfelberg
# calculate diffusion coeff in pore from traj file
# Units :
#   time : input ps ;
#   position : input A ;
#   D_out : A^2/ps

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pylab import *
import numpy as np
import os

wkDir = '/Users/lfelberg/PBSAM/pb_solvers/'\
                    'pbam/pbam_test_files/'\
                    'dynamics_test/contact_2sp/'
os.chdir(wkDir)
cgtot = 2

def trajDat( inFile):
    '''get data from traj file'''
    infile = open(inFile,"r")
    cgct = 0

    list_t, list_pos = [], [[] for x in range(cgtot)]
    DT, MSD, DR2 = 0, 0, 0
    for line in infile :
        line = line.split()
        if "Atom" in line[0]:
            dt = float(line[-1])
            list_t.append(dt)
        elif line[0] == "X" :
            pos_now = [np.array([ float(line[1]),
                                float(line[2]) , float(line[3]) ])]
            list_pos[cgct].append(pos_now)
            cgct += 1
            if cgct == cgtot:
                cgct = 0
        # ignore lines starting with 1

    list_t = np.array(list_t)
    list_pos = np.array(list_pos)

    return(list_t, list_pos)


fig = plt.figure()
ax = fig.add_subplot(111)
fontP = FontProperties()
fontP.set_size(18)

ct = 0
for dr in range(1,2, 1):
    dirN = '' #'0'+str(dr) +'/' if dr < 10 else str(dr) +'/'
    for traj in range(62):
        list_t, list_pos = trajDat(dirN+'dyn_contact_2sp_{0}.xyz'
                                    .format(traj))

        list_tuple, list_DT, list_DR2, list_DT2 = [], [], [], []
        list_DR22, list_varDT, list_varDR2 = [], [], []
        members = len(list_t)
        disp = np.zeros((cgtot, list_pos.shape[1]))

        for cg in range(cgtot):
            for step in range(list_pos.shape[1]):
                disp[cg][step] = np.linalg.norm(
                                list_pos[cg][0] - list_pos[cg][step])

        ax.plot(list_t/1000., disp[1], label=str(ct) ) # 'x', color='r', label="data")
        ct+=1
print ct

gcf().subplots_adjust(bottom=0.20, left=0.20)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(18)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)

ax.legend(prop=fontP)
ax.set_xlabel('Time (ns)', fontsize=24)
ax.set_ylabel('$\Delta r$', fontsize=24)
ax.set_title('')

plt.show()
plt.savefig("diffusion.pdf")

# calculate dt and mean square distance for slices of the data
# calculation of variance not working correctly yet
#   -- should have very small varDT and large varDR2
#for i in range(members, 1, -1) : # decrement by 1 until reaching 2
#    DT = 0
#    DT2 = 0
#    DR2 = 0
#    DR22 = 0
#    slices = members+1-i
#    for j in range(0, slices) :
#        arr_dt = list_t[0+j:i+j]
#        arr_pos = list_pos[0+j:i+j]
#        dt = arr_dt[-1] - arr_dt[0]
#        dr2 = np.sum( np.square( arr_pos[-1] - arr_pos[0] ) )
#        DT += dt
#        DT2 += dt*dt
#        DR2 += dr2
#        DR22 += dr2*dr2
#    list_DT.append(DT/slices)
#    list_DR2.append(DR2/slices)
#    if slices != 1 :
#        list_varDT.append( slices*(DT2 - DT*DT)/(slices*(slices-1)) )
#        list_varDR2.append( slices*( DR22 - DR2*DR2 )/(slices*(slices-1)) )
#    elif slices == 1 :
#        list_varDT.append(0)
#        list_varDR2.append(0)


# fit data for diffusion coeff calc.
# 1D or 3D ??? assuming 3D for now
# using dt from 3 to 15 ns
#   -- should make this depend on the data itself
#dim = 3
#
#list_y = []
## add directory for i/o
## reverse order s.t. smallest dt is first
#list_DT = list_DT[::-1]
#list_DR2 = list_DR2[::-1]
#
#fit = np.polyfit(list_DT[1:6], list_DR2[1:6], 1, full=True)
#D_out = fit[0][0] / ( 2 * dim )     # units of A^2/ps
#print(D_out)
#print(fit)
#for x in list_DT :
#    y = fit[0][0] * x + fit[0][1]
#    list_y.append(y)
#
## printing as text
#col1 = list_DT
#col2 = list_varDT
#col3 = list_DR2
#col4 = list_varDR2
#col5 = list_y
#output = '\n'.join('\t'.join(map(str,row)) for row in zip(col1,col3,col5))
#with open('diffusion.txt', 'w') as f:
#    f.write(output)

# for plots

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

'''
Program to plot a 3D version of the ESP from PB-AM
'''
dirName='/Users/felb315/'\
                 + 'Desktop/electrostatic_test/electro_barnase_test/0.05M/'
fileName = dirName + 'bars_0.05M_map.out'
#fileName = dirName + 'data/2fgr/2fgr_tri_move_map.out'
outFile= '/Users/felb315/Desktop/barnase_out.surf'
outFile= dirName + 'barnase_0.05M_3d_'


def FileOpen(fileName):
    """Gets data from 3D plot output of PB-AM"""
    lines = open(fileName).readlines()

    grid,org,dl = np.zeros(3), np.zeros(3), np.zeros(3)
    units, ct = 'jmol', 0
    pot = np.zeros((len(lines)-5, 4))

    for line in lines:
        temp = line.split()
        if 'units' in line[-10:-1]:
            units = temp[-1]
        elif 'grid' in line[0:10]:
            grid[0], grid[1] = int(temp[1]), int(temp[2])
        elif 'origin' in line[0:10]:
            org[0], org[1] = float(temp[1]), float(temp[2])
        elif 'delta' in line[0:10]:
            dl[0], dl[1] = float(temp[1]), float(temp[2])
        elif '#' not in line:
            temp = [float(x) for x in line.split()]

            pot[ct][0] = temp[0]
            pot[ct][1] = temp[1]
            pot[ct][2] = temp[2]
            pot[ct][3] = temp[3]
            ct += 1

    return(pot, org, dl, units)

def dispPlot( org, bn, xv, yv, zv, potential,
                title = '', lege = '', outFile = None ):
    """Plots the colormap of potential plot, 3D"""
    fig = plt.figure(1, figsize = (4, 4));
    ax = fig.add_subplot(111,projection='3d')

    n = len(xv)
    for i in range(n):
        if xv[i] < 0:
            break
    xv = xv[:i]
    yv = yv[:i]
    zv = zv[:i]
    potential=potential[:i]

    xv -= np.mean(xv)
    yv -= np.mean(yv)
    zv -= np.mean(zv)

    minl = min(min(xv), min(yv), min(zv))
    maxl = max(max(xv), max(yv), max(zv))

    big = max( abs(potential))
    cm = plt.get_cmap('seismic_r')
    cNorm = matplotlib.colors.Normalize(vmin=-big,
                                                               vmax=big)
    scalarMap = cmx.ScalarMappable(norm=cNorm,
                                                           cmap=cm)

    ax.scatter(xv, yv, zv, s = 85,
                     c=scalarMap.to_rgba(potential),
                     lw = 0)
    scalarMap.set_array(potential)
    #fig.colorbar(scalarMap)
    cbaxes = fig.add_axes([0.86, 0.1, 0.025, 0.75])
    cb = plt.colorbar(scalarMap, cax = cbaxes)
    cbaxes.set_xlabel(units, fontname='Arial',
            fontsize='small', labelpad=5.)

    ax.set_xlim([minl-2, maxl+2])
    ax.set_ylim([minl-2, maxl+2])
    ax.set_zlim([minl-2, maxl+2])

    #plt.title(title, fontsize = 13);
    ax.set_xlabel(r'$X (\AA)$', fontsize = 10)
    ax.set_ylabel(r'$Y (\AA)$', fontsize = 10)
    ax.set_zlabel(r'$Z (\AA)$', fontsize = 10 )
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    for tick in ax.zaxis.get_major_ticks():
        tick.label.set_fontsize(10)
    ax.tick_params(pad=-5.5)

    if outFile != None:
        for ii in range(0,360,180):
            ax.view_init(elev=20., azim=ii)
            #plt.gcf().subplots_adjust(bottom=0.15)
            plt.savefig(outFile+str(ii)+'.png',
                             bbox_inches='tight') #,dpi = 300)

    # plt.close()
    plt.show()

#------------------------------------------------------------------------------
# main

plt.close()
esp, org, dl, units = FileOpen(fileName)
if units == "jmol":
    units = "$J/$mol"
titl = 'Potential at surfaces in %s' % units

print(outFile)

dispPlot( org, dl, esp[:,0], esp[:,1], esp[:,2], esp[:,3],
          titl, outFile=outFile)

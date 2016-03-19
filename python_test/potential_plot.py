import numpy as np
import matplotlib.pyplot as plt

'''
Program to plot a 2D version of the ESP from PB-AM
'''

fileName='pot_x_1.00.dat'
outFile='pot_x_1.00.jpg'

#-----------------------------------------------------------------------

def FileOpen(fileName):
    lines = open(fileName).readlines()
    
    grid,org,dl = np.zeros(2), np.zeros(2),np.zeros(2)
    axval, ax, units = 0.0, 'x', 'jmol'
    pot = np.zeros((100,100))
    ct = 0

    for line in lines:
        temp = line.split()
        if 'units' in line[0:10]:
            units = temp[1]
        elif 'grid' in line[0:10]:
            grid[0], grid[1] = int(temp[1]), int(temp[2])  
            pot = np.zeros((grid[0], grid[1]))
        elif 'axis' in line[0:10]:
            ax, axval = temp[1], float(temp[2])
        elif 'origin' in line[0:10]:
            org[0], org[1] = float(temp[1]), float(temp[2])  
        elif 'delta' in line[0:10]:
            dl[0], dl[1] = float(temp[1]), float(temp[2])  
        elif '#' not in line:
            temp = [float(x) for x in line.split()]

            for i in range(int(grid[0])):
                pot[ct][i] = temp[i]
            ct += 1

    return(pot, org, dl, ax, axval, units)

def dispPlot( mn, bn, count, potential, title = '', 
                xlab = r'$X \AA$', ylab = r'$Y \, (\AA)$',
                lege = '', outFile = None ):
    fig = plt.figure(1, figsize = (4, 4)); 
    ax = fig.add_subplot(1,1,1)

    nbins = len(potential[0])
    #X, Y = np.meshgrid(np.arange(0, 100, 1), np.arange(0, 100, 1))
    #levels = np.arange(-1,1,0.1)
    #plt.contourf(X, Y, potential, levels)
    
    X = np.arange(mn[0], mn[0]+ nbins*bn[0], bn[0])
    Y = np.arange(mn[1], mn[1]+ nbins*bn[1], bn[1])
    plt.pcolor(X, Y, potential, vmin=-.10, vmax=.10)
    plt.colorbar()

    ax.set_xlim([X[0], X[-1]])
    ax.set_ylim([Y[0], Y[-1]])
    plt.title(title, fontsize = 13);
    ax.set_ylabel(ylab, fontsize = 10); ax.set_xlabel(xlab, fontsize = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    if outFile != None:
        plt.savefig(outFile,bbox_inches='tight', dpi = 200)
    plt.close()

#--------------------------------------------------------------------------------
# main

esp, org, dl, ax, axval, units = FileOpen(fileName)
boxl = len(esp[0])

titl = 'Cross section at ' + ax + ' = ' + str(axval)
titl += ' in ' + units
xla = r'$Y \, (\AA)$'; yla = r'$Z \, (\AA)$'
if ax == 'y':
    xla = r'$X \, (\AA)$'; yla = r'$Z \, (\AA)$'
elif ax == 'z':
    xla = r'$X \, (\AA)$'; yla = r'$Y \, (\AA)$'

dispPlot( org, dl, len(esp[0]), esp, titl,
                xlab = xla, ylab = yla,
                outFile=outFile)

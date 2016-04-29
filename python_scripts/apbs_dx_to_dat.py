import numpy as np
import math
import matplotlib.pyplot as plt

'''
Program to convert DX to a 2D version of the ESP 
'''
dirName= '/Users/lfelberg'
dirName += '/PBSAM/pb_solvers/python_test/APBS_compare/'
fileName= dirName + '1BRS.dx'
outFile=dirName + 'barn_apbs_z_0.00.jpg'
ind = 2 # chose 0 = x, 1 = y and 2 = z

barnCOM = [0, 0, 0]
barnRad2 = pow(25.6102,2)

#-----------------------------------------------------------------------
def FileOpen(fileName):
    """Gets data from DX output of APBS"""
    lines = open(fileName).readlines()
    
    grid,org,dl = np.zeros(3), np.zeros(3), np.zeros(3)
    pot = np.zeros((100,100, 100))
    xct, yct, zct = 0, 0, 0
    xpos, ypos, zpos = 0.0, 0.0, 0.0
    mx, mn, dlct = 0, 0, 0

    for line in lines:
        temp = line.split()
        if 'object 1' in line[0:10]:
            grid[0], grid[1] = int(temp[5]), int(temp[6]) 
            grid[2] = int(temp[7])
            pot = np.zeros((grid[0], grid[1], grid[2]))
            print np.shape(pot)
        elif 'origin' in line[0:10]:
            org[0], org[1] = float(temp[1]), float(temp[2])
            org[2] = float(temp[3]) 
        elif 'delta' in line[0:10]:
            dl[dlct] = float(temp[dlct+1]) 
            dlct += 1
        elif temp[0][0].isdigit() or temp[0][0] == '-':
            temp = [float(x) for x in line.split()]
            for i in range(len(temp)):
                mx = max( temp[i], mx)
                mn = min( temp[i], mn)
               
                # Finding actual pos in space
                xpos = xct*dl[0] + org[0]
                ypos = yct*dl[1] + org[1]
                zpos = zct*dl[2] + org[2]
                if ( pow(xpos-barnCOM[0],2) + 
                        pow(ypos-barnCOM[1],2) +
                        pow(zpos-barnCOM[2],2) ) < barnRad2:
                        temp[i] = float('nan')
                
                pot[xct][yct][zct] = temp[i]
                
                zct += 1
                if zct == grid[2]:
                    zct = 0
                    yct += 1
                    if yct == grid[1]:
                        yct = 0
                        xct += 1
    
    return(pot, org, dl, mx, mn)

def dispPlot( org, bn, count, potential, 
                mx = 0.1, mn = -0.1, title = '', 
                xlab = r'$X \AA$', ylab = r'$Y \, (\AA)$',
                lege = '', outFile = None ):
    """Plots the colormap of potential plot, 2D"""
    fig = plt.figure(1, figsize = (5, 4)); 
    ax = fig.add_subplot(1,1,1)

    nb = [len(potential[1]), len(potential)]
    
    X = np.arange(org[0], org[0]+ nb[0]*bn[0], bn[0])
    Y = np.arange(org[1], org[1]+ nb[1]*bn[1], bn[1])
    plt.pcolor(X, Y, potential, vmin=mn+0.1*mn, vmax=mx)
    plt.colorbar()

    ax.set_xlim([org[0], org[0]+ (nb[0]-1)*bn[0]])
    ax.set_ylim([org[1], org[1]+ (nb[1]-1)*bn[1]])
    plt.title(title, fontsize = 13);
    ax.set_ylabel(ylab, fontsize = 10)
    ax.set_xlabel(xlab, fontsize = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    if outFile != None:
        plt.savefig(outFile,bbox_inches='tight', dpi = 200)
    plt.close()
    
    return X, Y

#--------------------------------------------------------------------------------
# main

esp3D, org, dl, mx, mn  = FileOpen(fileName)

pos = barnCOM[ind]
idx = round((pos-org[ind]) / dl[ind])

ax = 'x'
xla = r'$Y \, (\AA)$'; yla = r'$Z \, (\AA)$'
esp = esp3D[idx, :, :]
mx, mn = 1.35, 0.0
if ind == 1:
    ax = 'y'
    esp = esp3D[:, idx, :]
    xla = r'$X \, (\AA)$'; yla = r'$Z \, (\AA)$'
    mx = 3.52078
    mx =  -0.393691
elif ind == 2:
    ax = 'z'
    esp = esp3D[:, :, idx]
    xla = r'$X \, (\AA)$'; yla = r'$Y \, (\AA)$'
    mx = 3.52078
    mx =  -0.393691

titl = "Cross section at {0} = {1:.2f}".format( ax,
                        org[ind] + dl[ind] * idx - barnCOM[ind])
titl += ' in kT'

org2d, dl2d = [], []
for i in range(3):
    if i != ind:
        org2d.append(org[i])    
        dl2d.append(dl[i])

X, Y = dispPlot( org2d, dl2d, len(esp[0]), np.transpose(esp), 
                mx , mn, title = titl, 
                xlab = xla, ylab = yla,
                outFile=outFile)
                


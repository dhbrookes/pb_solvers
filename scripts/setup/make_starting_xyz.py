import sys
import numpy as np
from numpy import linalg
import random

rad_avg = 27.0
nmol = [ ]

rand_rad = False
samp_grid = False

if ( sys.argv[1] == "sample_grid"):
    samp_grid = True
    rad_avg = float(sys.argv[3])
    ntypes = int(sys.argv[5])
    for typ in range(ntypes):
         nmol.append(int(sys.argv[6+typ]))
elif ( sys.argv[1] == "random_rad"):
    rand_rad = True
    radius = float(sys.argv[3])
else:
    print("Option not recognized!")

ntraj = int(sys.argv[2])
fname = sys.argv[4]

def sample_grid(rad_avg, nmol, nfiles, fname):
    ''' Sample a cubic grid of molecules for MFPT sim'''
    diam = 2.0*rad_avg;  pad = rad_avg/2.0
    ntypes = len(nmol); mol_tot = sum(nmol)

    bins = 0
    while bins**3 < mol_tot:
        bins += 1

    pos = np.zeros(( bins*bins*bins, 3))
    shift = (float(bins)/2.0) * (diam + pad)
    print bins, shift

    for nfile in range(nfiles):
        ct = 0
        for i in range(bins):
            for j in range(bins):
                for k in range(bins):
                    pos[ct] = np.array([i*diam+(i+1)*pad - shift,
                                                j*diam+(j+1)*pad - shift,
                                                k*diam+(k+1)*pad -shift])
                    ct += 1

        np.random.shuffle(np.array(pos))
        ct = 0
        for i in range(ntypes):
            f = open('{0}_{1:d}_{2:d}.xyz'.format(fname,i,nfile), 'w')
            for j in range(nmol[i]):
                f.write("{0:>8.2f}  {1:>8.2f}  {2:>8.2f}\n".format(
                        pos[ct][0], pos[ct][1], pos[ct][2]))

                ct += 1
            f.close()


def sample_spherical(npoints, ndim=3):
    '''Compute points on surface
    of a uniform sphere'''
    vec = np.random.randn(ndim, npoints)
    vec /= np.apply_along_axis(linalg.norm, 0, vec)
    return vec

def print_rand_to_xyz(sp_rad, nfile, xyz_pref):
    '''Print xyz files for points on sphere rad sp_rad'''
    xi, yi, zi = sample_spherical(nfile)

    for i in range(len(xi)):
        filen = xyz_pref + "_1_"+str(i+1) + ".xyz"
        f = open(filen, 'w')
        f.write("{0:>5.2f}  {1:>5.2f}  {2:>5.2f}".format(
                        sp_rad*xi[i],
                        sp_rad*yi[i],
                        sp_rad*zi[i]))
    f.close()

if (samp_grid):
    sample_grid(rad_avg, nmol, ntraj)

if (rand_rad):
    print_rand_to_xyz(radius, ntraj, fname)

import numpy as np
from numpy import linalg

ntypes = 2
rad_avg = 35.0
diam = 2.0*rad_avg
pad = rad_avg/2.0
nmol = [ 4, 4]
mol_tot = sum(nmol)

bins = 0
while bins**3 < mol_tot:
    bins += 1
    
pos = np.zeros(( bins*bins*bins, 3))

shift = (float(bins-1)/2.0) * (rad_avg + pad)

ct = 0
for i in range(bins):
    for j in range(bins):
        for k in range(bins):
            pos[ct] = np.array([i*diam+i*pad - shift,
                                          j*diam+j*pad - shift,
                                          k*diam+k*pad -shift])
            ct += 1

opts = np.arange(0, len(pos), 1)
for i in range(ntypes):
    samp = np.random.choice( opts, 
                                                size = nmol[i], 
                                                replace = False)                                   
    for j in range(nmol[i]):
        print pos[samp[j]][0], 
        print pos[samp[j]][1], pos[samp[j]][2]
        
    print ""
    np.delete(pos, samp)


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
        filen = xyz_pref + str(i+1) + ".xyz"
        f = open(filen, 'w')
        f.write("{0:5.3f} {1:5.3f} {2:5.3f}".format( 
                        sp_rad*xi[i], 
                        sp_rad*yi[i], 
                        sp_rad*zi[i]))
    f.close()


rand_rad = True

if (rand_rad):
    rad = 80.0
    print_rand_to_xyz(rad, 8, "pos_2_")


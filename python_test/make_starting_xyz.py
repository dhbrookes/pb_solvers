import numpy as np


ntypes = 2
rad_avg = 35.0
nmol = [ 4, 4]
mol_tot = sum(nmol)

diam = 2.0*rad_avg
pad = rad_avg/2.0

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
    
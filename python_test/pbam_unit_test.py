import numpy as np
import scipy 
from scipy.special import *
from scipy.misc import factorial

ns = np.arange(0.0,10.0,1.0)
zs = np.arange(1.0,11.0,9.0)

zs = np.array([0.16972088, 0.31519592])

resultsK = np.zeros((len(zs),len(ns)))
resultsI = np.zeros((len(zs),len(ns)))
nCt, zCt = 0, 0


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## for calculating the \hat{k}_n (z) of eq 3a, lotan 2006
for n in ns:
    zCt = 0
    for z in zs:
        mbfK = scipy.special.kv(n+0.5,z)
        knum = np.exp(z)*pow(z, n+0.5)
        kden = 1.0
        if n > 0:
            kden = float(reduce(int.__mul__,range(2*int(n)-1,0,-2)))
        resultsK[zCt][nCt] =  np.sqrt(2.0/np.pi)*(knum/kden)*mbfK
        
        zCt += 1
    nCt += 1
    
nCt, zCt = 0, 0
## for calculating the \hat{i}_n (z) of eq 3b, lotan 2006
for n in ns:
    zCt = 0
    for z in zs:
        mbfI = scipy.special.iv(n+0.5,z)
        inum = 1.0
        if n > 0:
            inum = float(reduce(int.__mul__,range(2*int(n)+1,0,-2)))
        iden = pow(z, n+0.5)
        resultsI[zCt][nCt] = np.sqrt(np.pi/2.0)*(inum/iden)*mbfI
        
        zCt += 1
    nCt += 1

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Useful SH and legendre constants:
nmax = 9.0
ms = np.arange(0.0,nmax,1.0)
const1 = ((2.0*nmax-1.0)/(nmax-ms))
const2 = ((ms+nmax-1.0)/(nmax-ms))

ls = np.arange(0.0,10.0,1.0)
dub_fac = np.zeros(len(ls))

for l in ls:
    dub_fac[l] = 1.0
    if l > 0:
        dub_fac[l] = reduce(int.__mul__,range(2*int(l)-1,0,-2))
        
## for calculating \sqrt{\frac{(n-|m|)!}{(n+|m|)!}}
nmax = 9.0
ms = np.arange(0.0,nmax+1.0,1.0)
sh_const = np.sqrt( factorial(nmax-ms) /
                                factorial(nmax+ms) )

## For calculating the associated Legendre Polynomial
## P_{n,m} (x) in eq 1 of Lotan 2006
theta = np.arange(0.0, np.pi + 0.1, np.pi/3.0)
n = 9.0

for z in theta:
    zCt = 0
    Pnm, deriv = scipy.special.lpmn(n, n, np.cos(z))
    #print z, np.cos(z)
    zCt += 1


## For calculating the associated spherical harmonics
## Y_{n,m} (\theta, \phi) in eq 1 of Lotan 2006
## In the paper they are different by a factor of
## \sqrt{ \frac{2n+1}{4\pi}}
theta = np.arange(0.0, np.pi + 0.1, np.pi/3.0)    # polar
phi = np.arange(0.0, 2.0*np.pi + 0.1, 0.5) # azimuthal
nmax = 10

theta = 0.5
phi = [0.5]

for p in phi:
    #print p
    for m in range(nmax):
        Ynm = scipy.special.sph_harm(m, nmax-1, p, theta)

       # print Ynm*pow(-1.0, m)*np.sqrt((4.0*np.pi) 
       #                                 /(2.0*float(nmax-1)+1.0))
        



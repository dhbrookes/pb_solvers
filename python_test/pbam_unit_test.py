import numpy as np
import scipy 
from scipy.special import *
from scipy.misc import factorial


bessel    = False
SHCons = False
SHCalc  = False
Rot         = True
nCt, zCt = 0, 0


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Calculating modified spherical bessel functions
if (bessel):
    ns = np.arange(0.0,10.0,1.0)
    zs = np.arange(1.0,11.0,9.0)
    #zs = np.array([0.16972088, 0.31519592])
    
    resultsK = np.zeros((len(zs),len(ns)))
    resultsI = np.zeros((len(zs),len(ns)))
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

if (SHCons):
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## For calculating the associated spherical harmonics
## Y_{n,m} (\theta, \phi) in eq 1 of Lotan 2006
## In the paper they are different by a factor of
## \sqrt{ \frac{2n+1}{4\pi}}

if (SHCalc or Rot):
    theta = np.arange(0.0, np.pi + 0.1, np.pi/3.0)    # polar
    phi = np.arange(0.0, 2.0*np.pi + 0.1, 0.5) # azimuthal
    nmax = 10
    
    for p in phi:
        for m in range(nmax):
            Ynm = scipy.special.sph_harm(m, nmax-1, p, theta)
    
            if SHCalc:
              print np.imag(Ynm*pow(-1.0, m)*np.sqrt((4.0*np.pi) 
                                          /(2.0*float(nmax-1)+1.0)))
        


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## For calculating rotation coefficients
## R_n^{0,s} = Y_{n,-s}(\theta)(\phi)

if (Rot):
    #theta = 0.0
    theta = 1.5953910356207
    #phi = [ 0.0 ]
    phi = 5.7258897721
    nmax = 5
    
    R = np.zeros(( 2*nmax, 2*nmax, 4*nmax),dtype=complex)
    
    for n in range(2*nmax):
        Ynm = scipy.special.sph_harm(range(0,n+1), n, phi, theta)
        for s in range(-n,n+1):
            if (s<0):
                val = Ynm[s+n]
            else:
                val = np.conj(Ynm[s])
            R[n][0][s+n] = val

    #print np.imag(Ynm*pow(-1.0, m)*np.sqrt((4.0*np.pi) 
    #                                    /(2.0*float(nmax-1)+1.0)))
                                        
    
    for m in range(0,nmax-1):
        for n in range(m+2,2*nmax-m):
            for s in range(-n+1, n):
                bmn = np.sqrt(float((n - m -1)*(n-m))/
                                        float((2.0*n-1.0)*2.0*n+1))
                
                bs1n = np.sqrt(float((n - ( s-1) -1)*(n-( s-1)))/
                                        float((2.0*n-1.0)*2.0*n+1))
                bs2n = np.sqrt(float((n - (-s-1) -1)*(n-(-s-1)))/
                                        float((2.0*n-1.0)*2.0*n+1))
                
                asn = np.sqrt(float((n+s+1)*(n-s+1))/float((2.0*n+1.0)*(2.0*n+3.0)))
                
                R[n][m+1][s+n] = (0.5*np.exp(complex(0.0,-1.0)*phi)*
                                                (1.0+np.cos(theta))*
                                                bs1n*R[n][m][s+n-1])
                R[n][m+1][s+n]-= (0.5*np.exp(complex(0.0, 1.0)*phi)*
                                            (1.0-np.cos(theta))*
                                            bs2n*R[n][m][s+n+1])
                R[n][m+1][s+n]+= np.sin(theta)*asn*R[n][m][s+n]
                R[n][m+1][s+n] *= (1.0/bmn)
            
    print R[4]

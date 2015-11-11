import numpy as np
import scipy 
from scipy.special import *
from scipy.misc import factorial

bessel    = False
SHCons = False
SHCalc  = False
Rot         = False
Trans      = True
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
        zCt += 1

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## For calculating the associated spherical harmonics
## Y_{n,m} (\theta, \phi) in eq 1 of Lotan 2006
## In the paper they are different by a factor of
## \sqrt{ \frac{2n+1}{4\pi}}

if (SHCalc):
    theta = np.arange(0.0, np.pi + 0.1, np.pi/3.0)    # polar
    phi = np.arange(0.0, 2.0*np.pi + 0.1, 0.5) # azimuthal
    nmax = 10
    
    for p in phi:
        for m in range(nmax):
            Ynm = scipy.special.sph_harm(m, nmax-1, p, theta)
    
        print Ynm*pow(-1.0, m)*np.sqrt((4.0*np.pi) 
                                          /(2.0*float(nmax-1)+1.0))
        


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## For calculating rotation coefficients
## R_n^{m,s} has a recursion relation
## R_n^{0,s} = Y_{n,-s}(\theta)(\phi)

if (Rot):
    theta = 0.0 
    phi =  0.0 
    #theta = 1.5953910356207
    #phi = 5.7258897721
    nmax = 10
    
    R = np.zeros(( 2*nmax, 2*nmax, 4*nmax),dtype=complex)
    
    A = np.zeros(( 2*nmax, 4*nmax))
    B = np.zeros(( 2*nmax, 4*nmax))
    
    for n in range(2*nmax):
        for m in range(-n, n+1):
            nD, mD = float( n ), float( m )
            A[n][m+nmax] =  np.sqrt(((nD+mD+1.)*(nD-mD+1.))
                                            /((2.*n+1.0)*(2.0*nD+3.)))
            sign = 1.0
            if (m<0):
                sign = -1.0
            B[n][m+nmax] =  sign * np.sqrt(((nD - mD -1)*(nD-mD))/
                                        ((2.*nD-1.)*(2.*nD+1.)))
                                        
    #for n in range(nmax):
    #    for m in range(-n, n+1):
    #        print( " {:.8f}, " .format( B[n][m+nmax])),     
    #    print ""
            
    for n in range(2*nmax):
        Ynm = scipy.special.sph_harm(range(0,n+1), n, phi, theta)
        for s in range(-n,n+1):
            if (s<0):
                val = Ynm[-s]
            else:
                val = np.conj(Ynm[s])
            R[n][0][s+2*nmax] = val*pow(-1.0, s)*np.sqrt((4.0*np.pi) 
                                          /(2.0*float(n)+1.0))
                                        
    
    for m in range(0,nmax):
        for n in range(m+2,2*nmax-m):
            for s in range(-n+1, n):
                
                sign1, sign2, sign3 = 1.0, 1.0, 1.0
                
                if (m<0):
                    sign1 = -1.0
                bmn = sign1 * np.sqrt(float((n - m -1)*(n-m))/
                                        float((2.0*n-1.0)*(2.0*n+1)))
                if ((s-1)<0):
                    sign2 = -1.0
                bs1n =  sign2 * np.sqrt(float((n - ( s-1) -1)*(n-( s-1)))/
                                        float((2.0*n-1.0)*(2.0*n+1)))
                                        
                if ((-s-1)<0):
                    sign3 = -1.0
                bs2n =  sign3 * np.sqrt(float((n - (-s-1) -1)*(n-(-s-1)))/
                                        float((2.0*n-1.0)*(2.0*n+1)))
                
                asn = np.sqrt(float((n+s+1)*(n-s+1))/float((2.0*n+1.0)*(2.0*n+3.0)))
                
                R[n-1][m+1][s+2*nmax] = (0.5*np.exp(complex(0.0,-1.0)*phi)*
                                                (1.0+np.cos(theta))*
                                                bs1n*R[n][m][s-1+2*nmax])

                R[n-1][m+1][s+2*nmax]-= (0.5*np.exp(complex(0.0, 1.0)*phi)*
                                            (1.0-np.cos(theta))*
                                            bs2n*R[n][m][s+1+2*nmax])
                R[n-1][m+1][s+2*nmax]+= np.sin(theta)*asn*R[n][m][s+2*nmax]
                R[n-1][m+1][s+2*nmax] *= (1.0/bmn)
    
    #n = nmax-1
    #for m in range(0,nmax): 
    #    print n ,m
    #    for s in range(0, n+1):
    #        print( "{:.6f}, " .format( float(np.real(R[n][m][s+2*nmax])))),
    #    
    #    print ("\n"),



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## For calculating rotation coefficients
## S_{n,l}^{m} has a recursion relation

if (Trans):
    #lam = 5.0
    z = 1.0
    #z = 8.132650244539
    lam = 5.0 #25.0
    nmax = 10
    
    kap = 0.0303073 #0.5
    
    S = np.zeros(( 2*nmax, 2*nmax, 4*nmax))
    
    alpha = np.zeros(( 2*nmax, 4*nmax))
    beta   = np.zeros(( 2*nmax, 4*nmax))
    nu      = np.zeros(( 2*nmax, 4*nmax))
    mu     = np.zeros(( 2*nmax, 4*nmax))
    
    for n in range(2*nmax):
        for m in range(-n, n+1):
            nD, mD = float( n ), float( m )
            
            alpha[n][m+nmax] = np.sqrt( ( nD + mD + 1.) * ( nD - mD + 1. ))
            beta[n][m+nmax]   = ( alpha[n][m+nmax] * pow( lam*kap, 2 ) ) / ( (2.*nD+1.)*(2.*nD+3.))
            
            sign = 1.0
            if (m<0):
                sign = -1.0
            nu[n][m+nmax] =  sign * np.sqrt(((nD - mD -1)*(nD-mD)))
            mu[n][m+nmax] = ( nu[n][m+nmax] * pow( lam*kap, 2 ) ) / ( (2*nD-1.)*(2.*nD+1.))
                                        
    #for n in range(nmax):
    #    for m in range(-n, n+1):
    #        print( " {:.8f}, " .format( mu[n][m+nmax])),     
    #    print ""
    
    
    resultsK = np.zeros(2*nmax)
    nCt = 0
    ## for calculating the \hat{k}_n (z) of eq 3a, lotan 2006
    for n in range(2*nmax):
        mbfK = scipy.special.kv(n+0.5,z*kap)
        knum = np.exp(z*kap)*pow(z*kap, n+0.5)
        kden = 1.0
        if n > 0:
            dfRan = range(2*n-1,0,-2)
            for i in dfRan:
                kden *= i
        resultsK[n] =  np.sqrt(2.0/np.pi)*(knum/kden)*mbfK
        
        S[0][n][0] = pow(lam/z, n) * np.exp(-kap*z) * (1.0/z) * resultsK[n] 
        S[n][0][0] = pow(-1.0, n) * S[0][n][0]
        
    for n in range(nmax): 
        print( "{:.6f}, ".format(S[n][0][0])), 
    print ""
    
    for l in range(1, 2*nmax - 2 ):
        a00  = np.sqrt( float( ( (0) + 0 + 1) * ( (0) - 0 + 1 )))
        al0   = np.sqrt( float( ( ( l) + 0 + 1) * ( ( l) - 0 + 1 )))  
        al10  = np.sqrt( float( ( (l-1) + 0 + 1) * ( (l-1) - 0 + 1 )))     
        bl10 = ( al10 * pow( lam*kap, 2 ) ) / ( float( (2*(l-1)+1)*(2*(l-1)+3)))    
        
        S[1][l][0] = (-1.0/a00)*(bl10*S[0][l-1][0] + al0*S[0][l+1][0])
    
    for n in range(1, nmax-1):
        for l in range(n+1, 2*nmax - n - 2 ):
            al0 = np.sqrt( float( ( (l-1) + 0 + 1) * ( (l-1) - 0 + 1 )))
            bl0 = ( al0 * pow( lam*kap, 2 ) ) / ( float( (2*(l-1)+1)*(2*(l-1)+3))) 
            
            an0 = np.sqrt( float( ( (n-1) + 0 + 1) * ( (n-1) - 0 + 1 )))
            bn0 = ( an0 * pow( lam*kap, 2 ) ) / ( float( (2*(n-1)+1)*(2*(n-1)+3))) 
            
            al = np.sqrt( float( ( (l) + 0 + 1) * ( (l) - 0 + 1 )))
            an = np.sqrt( float( ( (n) + 0 + 1) * ( (n) - 0 + 1 )))
            
            S[n+1][l][0]   = bl0 * S[n][l-1][0] + bn0 * S[n-1][l][0]
            S[n+1][l][0] += al * S[n][l+1][0]
            S[n+1][l][0] *=  ( -1.0 / an ) 
    
    #for n in range(nmax): 
    #    print( "{:.6f}, ".format(S[n][n+1][0])),  

    for m in range(1, nmax):
        for l in range(m, 2*nmax - m - 1):
            
            sign1, sign2, sign3 = 1.0, 1.0, 1.0
            
            if ( (-m) < 0 ):
                sign1 = -1.0
                sign2 = -1.0
                
            if ( (m-1) < 0 ):
                sign3 = -1.0
                
            eta1 = sign1 * np.sqrt( float( ( m - (-m) - 1.0 ) * (m - (-m) )))
            
            eta2 = sign2 * np.sqrt( float( ( l - (-m) - 1.0 ) * (l - (-m) )))
            mu2 = ( eta2 * pow( lam*kap, 2 ) ) / float( (2*l-1)*(2*l+1))
            
            eta3 = sign3 * np.sqrt( float( ( (l+1) - (m-1) - 1.0 ) * ((l+1) - (m-1) )))
            
            S[m][l][m] = mu2 * S[m-1][l-1][m-1] + eta3 * S[m-1][l+1][m-1]
            S[m][l][m] *=  ( -1.0 / eta1 ) 
            
    #for m in range(0, nmax): 
    #    print( "{:.8f}, ".format(S[m][m+2][m])),  
    
    for m in range(1, nmax):
        for n in range( m, nmax-1):
            for l in range( n+1, 2*nmax - n - 2):
                
                al1m = np.sqrt( float( ( (l-1) + m + 1) * ( (l-1) - m + 1 )))
                bl1m = ( al1m * pow( lam*kap, 2 ) ) / ( float( (2*(l-1)+1)*(2*(l-1)+3))) 
            
                an1m = np.sqrt( float( ( (n-1) + m + 1) * ( (n-1) - m + 1 )))
                bn1m = ( an1m * pow( lam*kap, 2 ) ) / ( float( (2*(n-1)+1)*(2*(n-1)+3))) 
                
                alm = np.sqrt( float( ( l + m + 1) * ( l - m + 1 )))
                
                anm = np.sqrt( float( ( n + m + 1) * ( n - m + 1 )))
                
                S[n+1][l][m]  = bl1m * S[n][l-1][m] + bn1m * S[n-1][l][m] + alm* S[n][l+1][m]
                S[n+1][l][m] *= ( -1.0 / anm )
    #for m in range(1, nmax-1): 
    #    print( "{:.8f}, ".format(S[m+2][m+3][m])),  
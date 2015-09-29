//
//  SHCalc.hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef SHCalc_h
#define SHCalc_h

#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
#include "MyMatrix.h"

using namespace std;


/*
 Class for storing the constants that may be used in multiple spherical harmonics calculations
 */
class SHCalcConstants
{
protected:
    int                             nPoles_;  // number of poles
    MyMatrix<double>                legConsts1_;  // (2l-1)/(l-m) for use in legendre computation
    MyMatrix<double>                legConsts2_;  // (l+m-1)/(l-m) for use in legendre computation
    MyMatrix<double>                shConsts_;  // sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
    vector<double>                  dubFac_;  // (2l-1)!! double factorial, for use in legendre recursion
    
public:
    SHCalcConstants(const int N);
    
    const double get_leg_consts1_val(const int n, const int m) const    { return legConsts1_(n, m); }
    const double get_leg_consts2_val(const int n, const int m) const    { return legConsts2_(n, m); }
    const double get_sh_consts_val(const int n, const int m) const      { return shConsts_(n, m);   }
    const double get_dub_fac_val(const int i) const                     { return dubFac_[i];        }
    const int get_n() const                                             { return nPoles_;                }
};

/*
 Class for computing spherical harmonics. This includes
 calcualtion of the associated Legendre Polynomials
 
 The spherical harmonics in this case are defined by the equation:
 Y_(n,m)(theta, phi) = (-1)^m * sqrt((n-m)! / (n + m)!) * P_(n,m)(cos(theta)) * exp(i*m*phi)
 where P_(n, m) are the associated Legendre polynomials.
 
 These are constrcucted dynamically and returned as a matrix of values for every n,m
 */


class SHCalc
{
protected:
    
    int                         nPoles_;  //number of poles (the output matrix will be 2Nx2N)
    double                      theta_;
    double                      phi_;
    const SHCalcConstants*      _consts_;
    MyMatrix<double>            P_; // legendre polynomial calculated for this spherical harmonics
    MyMatrix<complex<double> >  Y_;  // the spherical harmonics calcualted by this class
    
    void calc_legendre();  // calculate the legendre polynomial at every n, m (store in this.P_)
    void calc_sh(); // calculate the spherical harmonics at every n, m  (store in this.Y_)

public:
    
    /*
     Constructor given size and constants:
     */
    SHCalc(const int N, const SHCalcConstants* _consts,
           const double theta, const double phi);
    
    
    const complex<double> get_result(const int n, const int m) const;  // retrieve the result for n, m values
    
};

#endif /* SHCalc_hpp */

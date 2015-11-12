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
#include "Constants.h"
#include "MyMatrix.h"
#include "util.h"

using namespace std;


/*
 Class for storing the constants that may be used in
 multiple spherical harmonics calculations
 */
class SHCalcConstants
{
protected:
  int                numVals_;    // number of poles
  MyMatrix<double>   legConsts1_; // (2l-1)/(l-m) for use in legdre comp
  MyMatrix<double>   legConsts2_; // (l+m-1)/(l-m) for use in legdre compn
  MyMatrix<double>   shConsts_;  // sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
  vector<double>     dubFac_;    // (2l-1)!! double factrl, in legendre recursion
  
public:
  SHCalcConstants(const int num_vals=Constants::MAX_NUM_POLES);
  
  const double get_leg_consts1_val(const int n, const int m)
  { return legConsts1_(n, m); }
  const double get_leg_consts2_val(const int n, const int m)
  { return legConsts2_(n, m); }
  const double get_sh_consts_val(const int n, const int m)
  { return shConsts_(n, m);   }
  const double get_dub_fac_val(const int i) const
  { return dubFac_[i];        }
  const int get_n() const
  { return numVals_;           }
};


/*
 Class for computing spherical harmonics. This includes
 calcualtion of the associated Legendre Polynomials
 
 The spherical harmonics in this case are defined by the equation:
 Y_(n,m)(theta, phi) = (-1)^m * sqrt((n-m)! / (n + m)!) *
 P_(n,m)(cos(theta)) * exp(i*m*phi)
 where P_(n, m) are the associated Legendre polynomials.
 
 These are constructed dynamically and returned as a matrix of
 values for every n,m
 */
class SHCalc
{
protected:
  
  int                     numVals_;  //# of poles (output matrix will be 2Nx2N)
  shared_ptr<SHCalcConstants> _consts_;
  MyMatrix<double>        P_;  // legendre polynomials
  MyMatrix<cmplx>         Y_;  // spherical harmonics calcualted by this class
  
  void calc_legendre(const double theta);  // calculate the legendre polynomial
  //  at every n, m (store in this.P_)
  
public:
  
  SHCalc(const int num_vals, SHCalcConstants consts);
  
  // calculate the spherical harmonics at every n, m  (store in this.Y_)
  void calc_sh(const double theta, const double phi);
  
  // retrieve the result for n, m values
  cmplx get_result(const int n, const int m);
  
  // retrieve the full calculated Y_ matrix
  MyMatrix<cmplx> get_full_result() { return Y_; }
  
  double get_legendre_result( int n, int m );
};

#endif /* SHCalc_hpp */
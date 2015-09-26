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
#include "MyMatrix.h"

using namespace std;

/*
 Class for computing spherical harmonics. This includes
 calcualtion of the associated Legendre Polynomials
 */


class SHCalc
{
protected:
    
    int                             N_;
    MyMatrix<double>                consts1_;  // (2l-1)/(l-m) for use in legendre computation
    MyMatrix<double>                consts2_;  // (l+m-1)/(l-m) for use in legendre computation
    MyMatrix<double>                consts3_;  // sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
    vector<double>                  consts4_;  // (2l-1)!! double factorial, for use in legendre recursion
    vector<int>                     IDX_;  //<! An index vector for 2*Pol + 1 values
//    double          RS_;  //!< Scaling factor for system (~ protein radius)
//    double          IRS_;  //!< Inverse scaling factor (=1/RS)
//    double          kappa_;  //!< Inverse Debye length
    
    vector<complex<double> > mM_;  //!< A vector of complex numbers of matrix coefficients
    vector<double> mP_;  //!< Vector of legendre polynomials of NPOLES * NPOLES
    
public:
    
    /*
     Initialize constants
     */
    SHCalc(int N);
    
    
    
};


#endif /* SHCalc_hpp */

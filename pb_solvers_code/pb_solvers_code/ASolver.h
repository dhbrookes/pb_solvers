//
//  ASolver.hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef ASolver_h
#define ASolver_h

#include <stdio.h>
#include "MyMatrix.h"
#include "BesselCalc.h"
#include "Constants.h"
#include "SHCalc.h"
#include "System.h"
#include "util.h"

/*
 This class is designed to compute the vector A defined on pg. 544 of Lotan, Head-Gordon 2006
 */
class ASolver
{
protected:
    typedef MyMatrix<MyMatrix<double> > MatOfMats;
    typedef MyVector<MyVector<double> > VecOfVecs;
    typedef MyVector<MyMatrix<double> > VecOfMats;
    
    int                   N_;  // number of molecules
    int                   p_;  // max value for n (2*nPoles_ usually)
    VecOfVecs             A_;
    VecOfMats             E_;
    MatOfMats             gamma_, delta_, T;
    const BesselCalc*     _besselCalc_;
    const SHCalc*         _shCalc_;
    const System          sys_;  // system data (radii, charges, etc.)
    const Constants       consts_;
    
    
    const double calc_indi_gamma(int i, int n) const;  // calculate one index of inner gamma matrix
    const double calc_indi_delta(int i, int n) const;  // calculate on index of inner delta matrix
    const double calc_indi_e(int i, int n, int m);
    
    void compute_gamma();  // compute the gamma matrix (as defined on page 544 of Lotan(2006)
    void compute_delta();  //comput the delta
    void compute_E();  // compute the E vector (equations of page 543 of Lotan 2006)

public:
    
    ASolver(const vector<double>* a, const int N, const int p,
            const BesselCalc* _bcalc, const SHCalc* _shCalc,
            const System sys);
    
};

#endif /* ASolver_hpp */

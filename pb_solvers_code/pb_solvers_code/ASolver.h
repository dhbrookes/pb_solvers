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


/*
 This class is designed to compute the vector A defined on pg. 544 of Lotan, Head-Gordon 2006
 */
class ASolver
{
protected:
    typedef MyMatrix<MyMatrix<double> > MatOfMats;
    typedef MyVector<MyVector<double> > VecOfVecs;
    
    int                     N_;  // number of molecules
    int                     p_;  // max value for n (2*nPoles_ usually)
    VecOfVecs               A_, E_;
    MatOfMats               gamma_, delta_, T;
    const BesselCalc*       _besselCalc_;
    const Constants*        _consts_;
    const vector<double>*   _a_; // vector of molecular radii
    
    
    const double calc_indi_gamma(int i, int n) const;  // calculate one index of inner gamma matrix
    const double calc_indi_delta(int i, int n) const;  // calculate on index of inner delta matrix
    
    void compute_gamma();  // compute the gamma matrix (as defined on page 544 of Lotan(2006)
    void compute_delta();  //comput the delta

public:
    
    ASolver(const vector<double>* a, const int N, const int p,
            const BesselCalc* _bcalc, const Constants* _consts);
    
};

#endif /* ASolver_hpp */

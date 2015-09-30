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
 This class is designed to compute the vector A defined in Equation 22 -- Lotan 2006, page 544
 */
class ASolver
{
protected:
    
    int                         N_;  // number of molecules
    int                         p_;  // max value for n (2*nPoles_ usually)
    VecOfMats<cmplx>::type      E_;
    MatOfMats<double>::type     gamma_, delta_;
    const BesselCalc*           _besselCalc_;
    const System                sys_;  // system data (radii, charges, etc.)
    const Constants             consts_;
    SHCalc*                     _shCalc_;
    
    // pre-computed spherical harmonics matrices for every charge in the system
    // inner vector is all SH for all the charges in a molecule. Outer vector is every molecule
    vector<vector<MyMatrix<cmplx> > > all_sh;
    
    vector<MyMatrix<cmplx> > calc_mol_sh(Molecule mol); // calculate the SH for all charges in a molecule
    const double calc_indi_gamma(int i, int n) const;  // calculate one index of inner gamma matrix
    const double calc_indi_delta(int i, int n) const;  // calculate on index of inner delta matrix
    const cmplx calc_indi_e(int i, int n, int m);
    
    void compute_gamma();  // compute the gamma matrix (as defined on page 544 of Lotan 2006)
    void compute_delta();  //comput the delta
    void compute_E();  // compute the E vector (equations on page 543 of Lotan 2006)

public:
    
    ASolver(const vector<double>* a, const int N, const int p,
            const BesselCalc* _bcalc, SHCalc* _shCalc,
            const System sys);
    
};

#endif /* ASolver_hpp */

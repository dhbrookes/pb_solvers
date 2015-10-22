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
#include "ReExpCalc.h"
#include "SHCalc.h"
#include "System.h"
#include "util.h"

/*
 This class is designed to compute the vector A defined in Equation 22 
 Lotan 2006, page 544
 */
class ASolver
{
protected:
    
  int                         N_;  // number of molecules
  int                         p_;  // max value for n (2*numVals_ usually)    
  VecOfMats<cmplx>::type      E_;
  MatOfMats<double>::type     gamma_, delta_;
  const BesselCalc*           _besselCalc_;
  System*                     _sys_;  // system data (radii, charges, etc.)
  const Constants*            _consts_;
  SHCalc*                     _shCalc_;
  ReExpCoeffsConstants*       _reExpConsts_;
  
//  MyMatrix<ShPt> interDists_; //inter molecular vector
//  // pre-computed spherical harmonics for every inter molecular vector
//  MatOfMats<cmplx>::type all_inter_sh;
  
  //re expansion coefficients calcualted for every inter molecular vector
  MyMatrix<ReExpCoeffs>  T_;
  
  // pre-computed spherical harmonics matrices for every charge in the system
  // inner vector is all SH for all the charges in a molecule.
  // Outer vector is every molecule
  vector<vector<MyMatrix<cmplx> > > all_sh;
  
  // calculate the SH for all charges in a molecule
  vector<MyMatrix<cmplx> > calc_mol_sh(Molecule mol);
  
  // calculate one index of inner gamma matrix
  double calc_indi_gamma(int i, int n);
  
  // calculate on index of inner delta matrix
  double calc_indi_delta(int i, int n);
  cmplx calc_indi_e(int i, int n, int m);
  
  // compute the gamma matrix (as defined on page 544 of Lotan 2006):
  void compute_gamma();
  
  //compute the delta matrix (as defined on page 544 of Lotan 2006):
  void compute_delta();
  
  // compute the E vector (equations on page 543 of Lotan 2006)
  void compute_E();

public:
  
  MatOfMats<double>::type get_gamma()   { return gamma_; }
  MatOfMats<double>::type get_delta()   { return delta_; }
  VecOfMats<cmplx>::type  get_E()       { return E_; }
  
  double get_gamma_ni( int i, int n)    { return gamma_( i, i )( n, n); }
  double get_delta_ni( int i, int n)    { return delta_( i, i )( n, n); }
  cmplx  get_E_ni( int i, int n, int m) { return E_[ i ]( n, m+n ); }
  
  ASolver(const int N, const int p, const BesselCalc* _bcalc,
          SHCalc* _shCalc, System* sys);
  virtual ~ASolver();
  ASolver(const ASolver& other);
  ASolver& operator=(const ASolver& other);
  
}; // End ASolver

#endif /* ASolver_hpp */

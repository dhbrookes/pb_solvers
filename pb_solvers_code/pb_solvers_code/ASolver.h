//
//  ASolver.h
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
  
  VecOfMats<cmplx>::type      A_;  // what we iteratively solve for
  
  int                         N_;  // number of molecules
  int                         p_;  // max value for n (2*numVals_ usually)    
  VecOfMats<cmplx>::type      E_;
  VecOfMats<cmplx>::type      gamma_, delta_;
  
  ReExpCoeffsConstants        reExpConsts_;

  shared_ptr<BesselCalc>      _besselCalc_;
  shared_ptr<System>          _sys_;  // system data (radii, charges, etc.)
  shared_ptr<SHCalc>          _shCalc_;
  
  // re expansion coefficients calculated for every inter molecular vector
  MyMatrix<ReExpCoeffs>  T_;
  
  // pre-computed spherical harmonics matrices for every charge in the system
  // inner vector is all SH for all the charges in a molecule.
  // Outer vector is every molecule
  vector<vector<MyMatrix<cmplx> > > all_sh;
  
  // calculate the SH for all charges in a molecule
  vector<MyMatrix<cmplx> > calc_mol_sh(Molecule mol);
  
  // calculate one index of inner gamma matrix
  cmplx calc_indi_gamma(int i, int n);
  
  // calculate on index of inner delta matrix
  cmplx calc_indi_delta(int i, int n);
  cmplx calc_indi_e(int i, int n, int m);
  
  // pre-compute spherical harmonics matrices for every charge in the system
  void pre_compute_all_sh();
  
  // Compute the T matrix (re expansion coefficients for
  // every inter molecular vector)
  void compute_T();
  
  // compute the gamma matrix (as defined on page 544 of Lotan 2006):
  void compute_gamma();
  
  // compute the delta matrix (as defined on page 544 of Lotan 2006):
  void compute_delta();
  
  // compute the E vector (equations on page 543 of Lotan 2006)
  void compute_E();
  
  // initialize A vector
  void init_A();
  
  // re-expand element i of A withh element (i, j) of T and return results
  MyMatrix<cmplx> re_expandA(int i, int j);
  
  // perform one iteration of the solution for A (eq 51 in Lotan 2006)
  void iter();
  

public:
  
  VecOfMats<cmplx>::type&  get_gamma()   { return gamma_; }
  VecOfMats<cmplx>::type&  get_delta()   { return delta_; }
  VecOfMats<cmplx>::type&  get_E()       { return E_; }
  VecOfMats<cmplx>::type&  get_A()       { return A_; }
  
  cmplx get_gamma_ni( int i, int n)    { return gamma_[i]( n, n); }
  cmplx get_delta_ni( int i, int n)    { return delta_[i]( n, n); }
  cmplx get_E_ni( int i, int n, int m) { return E_[ i ]( n, m+p_ ); }
  cmplx get_A_ni(int i, int n, int m)  { return A_[ i ]( n, m+p_ ); }
  
  //void set_E_ni( int i, int n, int m)  { E_[ i ]( n, m+p_ ); }
  //void set_A_ni( int i, int n, int m)  { A_[ i ]( n, m+p_ ); }
  
  ASolver(const int N, const int p, BesselCalc bcalc,
          SHCalc shCalc, System sys);

  
  
  //numerically solve for A given the number of desired iterations
  void solve_A(int num_iter);
  
}; // End ASolver

#endif /* ASolver_h */

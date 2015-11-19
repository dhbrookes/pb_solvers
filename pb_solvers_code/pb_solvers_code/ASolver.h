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
#include <iomanip>
#include <iostream>
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
  
  VecOfMats<cmplx>::type      A_;  // solution
  VecOfMats<cmplx>::type      prevA_;  // previous value for
                                       // calculating convergence criteria
  
  bool solvedA_;
  
  /*
   Gradient of A. Each (i, j) entry in the outer vector is grad_j(A^(i))
   The inner vectors are all of length three and represent dA/dr, dA/dtheta
   and dA/dphi, respectively. This can only be calculated once A has been
   solved for
   */
  MyMatrix<VecOfMats<cmplx>::type > gradA_;
  MyMatrix<VecOfMats<cmplx>::type > prevGradA_;
  
  
  int                         N_;  // number of molecules
  int                         p_;  // max value for n (2*numVals_ usually)
  double                  a_avg_;  // the average radius of particles in syst
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
  
  // inialize grad(A) matrix to the zero matrix
  void init_gradA();
  
  /*
   enum for telling the ReExpCoeffs which values to retrieve 
   in the below methods
   */
  enum WhichReEx { BASE, DDR, DDPHI, DDTHETA };
  
  // re-expand element j of A with element (i, j) of T and return results
  // if prev=True then re-expand prevA
  MyMatrix<cmplx> re_expandA(int i, int j, bool prev=false);
  
  
  
  // re-expand element j of grad(A) with element (i, j) of T
  VecOfMats<cmplx>::type re_expand_gradA(int i, int j, int wrt, bool prev=false);
  
  // re-expand element j of A with element (i, j) of grad(T) and return results
  // uses eq 46 to solve eq 47 in Lotan 2006
  VecOfMats<cmplx>::type re_expandA_gradT(int i, int j, bool prev=false);
  
  // perform first part of T*A and return results (see eq 46 in Lotan 2006)
  // input wrt is only used if whichA is not BASE (then we need to know which
  // molecule the gradient is with respect to). If prev=True then expand
  // prevA (or prevGradA_)
  MyMatrix<cmplx> expand_RX(int i, int j, WhichReEx whichR,
                            WhichReEx whichA, bool prev, int wrt=-1);
  
  // perform second part of T*A and return results (see eq 46 in Lotan 2006)
  MyMatrix<cmplx> expand_SX(int i, int j, MyMatrix<cmplx> x1,
                            WhichReEx whichS);
  
  // perform third part of T*A and return results (see eq 46 in Lotan 2006)
  MyMatrix<cmplx> expand_RHX(int i, int j, MyMatrix<cmplx> x2,
                             WhichReEx whichRH);
  
  // convenience method for retrieving values from A and gradA (or their
  // previous values
  cmplx which_aval(WhichReEx whichA, bool prev, int i, int n,
                   int m, int wrt=-1);

  // perform one iteration of the solution for A (eq 51 in Lotan 2006)
  void iter();
  
  // perform one iterations of the solution for grad(A) (eq53 in Lotan 2006)
  void grad_iter();
  
  // calculate the change in A_ from prevA_ (eq 52 in Lotan 2006)
  double calc_change(WhichReEx whichA=BASE, int wrt=-1);
  
  // sum of many calls to the above
  double calc_grad_change();

public:
  
  VecOfMats<cmplx>::type&  get_gamma()     { return gamma_; }
  VecOfMats<cmplx>::type&  get_delta()     { return delta_; }
  VecOfMats<cmplx>::type&  get_E()         { return E_; }
  VecOfMats<cmplx>::type&  get_A()         { return A_; }
  
  VecOfMats<cmplx>::type calc_L();
  
  cmplx get_gamma_ni( int i, int n)        { return gamma_[i]( n, n); }
  cmplx get_delta_ni( int i, int n)        { return delta_[i]( n, n); }
  
  cmplx get_SH_ij(int i, int j, int n, int m)
                                           { return all_sh[i][j]( n, abs(m)); }
  cmplx get_E_ni(int i, int n, int m)      { return E_[ i ]( n, m+p_ ); }
  //get element of A^(i)
  cmplx get_A_ni(int i, int n, int m)      { return A_[ i ]( n, m+p_ ); }
  cmplx get_prevA_ni(int i, int n, int m)  { return prevA_[ i ]( n, m+p_ ); }
  
  
  // get elements of grad_j(A^(i))
  cmplx get_dAdr_ni(int i, int j, int n, int m)
  {
    return gradA_(i, j)[0](n, m+p_);
  }
  cmplx get_dAdtheta_ni(int i, int j, int n, int m)
  {
    return gradA_(i, j)[1](n, m+p_);
  }
  cmplx get_dAdphi_ni(int i, int j, int n, int m)
  {
    return gradA_(i, j)[2](n, m+p_);
  }
  
  cmplx get_prev_dAdr_ni(int i, int j, int n, int m)
  {
    return prevGradA_(i, j)[0](n, m+p_);
  }
  cmplx get_prev_dAdtheta_ni(int i, int j, int n, int m)
  {
    return prevGradA_(i, j)[1](n, m+p_);
  }
  cmplx get_prev_dAdphi_ni(int i, int j, int n, int m)
  {
    return prevGradA_(i, j)[2](n, m+p_);
  }
  
  void set_A_ni(int i, int n, int m, cmplx val)
  {
    (&A_[i])->set_val( n, m+p_, val);
  }
  
  ASolver(const int N, const int p, BesselCalc bcalc,
          SHCalc shCalc, System sys);

  void print_Ei( int i, int p);
  
  void print_Ai( int i, int p);


  //numerically solve for A given the desired precision
  void solve_A(double prec);
  
  // numerically solve for grad(A) given the desired precision
  // must solve for A before this
  void solve_gradA(double prec);
  
}; // End ASolver

#endif /* ASolver_h */

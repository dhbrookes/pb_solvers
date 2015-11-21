//
//  EnergyForce.h
//  pb_solvers_code
//
//  Created by David Brookes on 11/18/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef EnergyForce_h
#define EnergyForce_h

#include <stdio.h>
#include <memory>
#include "util.h"
#include "MyMatrix.h"
#include "ASolver.h"

using namespace std;

/*
 Class for calculating the energy of molecules in the system given 
 an ASolver object
 */
class EnergyCalc
{
protected:
  shared_ptr<VecOfMats<cmplx>::type> _A_;
  shared_ptr<VecOfMats<cmplx>::type> _L_;
  double epsS_;  // solvent dielectric constant
  int N_;  // number of molecules
  int p_;  // max number of poles
  
  MyVector<cmplx> omega_;  // result of energy calculation
  
public:
  
  EnergyCalc(VecOfMats<cmplx>::type A, VecOfMats<cmplx>::type L,
             double epsS, int N, int p);
  
  // fill omega_
  void calc_energy();
  
  
  // get the energy for a specific molecule:
  cmplx get_omega_i(int i)  { return omega_[i]; }
  // get all energy:
  MyVector<cmplx> get_omega() { return omega_; }
};

#endif /* EnergyForce_h */

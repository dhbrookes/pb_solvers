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
  Constants const_;
  
  /*
   Enum for the units of energy
   */
  enum WhichUnit { INTER, KCALMOL, KT, JMOL };
  
  MyVector<double> omega_;  // result of energy calculation, internal units
  
public:
  EnergyCalc(VecOfMats<cmplx>::type A, VecOfMats<cmplx>::type L,
             Constants const_, int N, int p);
  
  // fill omega_
  void calc_energy();
  
  // get the energy for a specific molecule:
  double get_omega_i_int(int i)  { return omega_[i]; }
  // get all energy:
  MyVector<double> get_omega_int() { return omega_; }
  
  // energy in kCal/mol:
  double get_omega_i_kcal(int i)
  { return const_.convert_int_to_kcal_mol(omega_[i]); }
  MyVector<double> get_omega_kcal()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = const_.convert_int_to_kcal_mol(omega_[n]);
    return omeg;
  }
  
  // energy in kT:
  double get_omega_i_kT(int i)
  { return const_.convert_int_to_kT(omega_[i]); }
  MyVector<double> get_omega_kT()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = const_.convert_int_to_kT(omega_[n]);
    return omeg;
  }
  
  // energy in joules/mol:
  double get_omega_i_jmol(int i)
  { return const_.convert_int_to_jmol(omega_[i]); }
  MyVector<double> get_omega_jmol()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = const_.convert_int_to_jmol(omega_[n]);
    return omeg;
  }
};

#endif /* EnergyForce_h */

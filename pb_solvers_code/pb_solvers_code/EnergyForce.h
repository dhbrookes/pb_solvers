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

/*
 Class for calculating the forces on molecules in the system given
 an ASolver object
 */
class ForceCalc
{
protected:
  shared_ptr<VecOfMats<cmplx>::type> _A_;
  shared_ptr<VecOfMats<cmplx>::type> _L_;
  
  shared_ptr< MyMatrix<VecOfMats<cmplx>::type > > _gradA_;
  shared_ptr< MyVector<VecOfMats<cmplx>::type > > _gradL_;
  
  double epsS_;
  int N_;
  int p_;
  Constants const_;
  
  VecOfVecs<double>::type F_;
  
public:
  ForceCalc(VecOfMats<cmplx>::type A,
            MyMatrix<VecOfMats<cmplx>::type > gradA_,
            VecOfMats<cmplx>::type L,
            MyVector<VecOfMats<cmplx>::type > gradL_,
            Constants con, int N, int p);
  
  void calc_force();  // fill F_
  
  MyVector<double> get_fi(int i)     { return F_[i]; }
  VecOfVecs<double>::type get_F()    { return F_; }
  
};

/*
 Class for calculating the torque on every molecule in the system
 */
class TorqueCalc
{
protected:
  
  // outer vector has an entry for every molecule. Inner vector is the torque
  // on that molecule
  VecOfVecs<cmplx>::type tau_;
  
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr< MyVector<VecOfMats<cmplx>::type > > _gradL_;
  
  Constants consts_;
  shared_ptr<System> _sys_;
  
  double epsS_;
  int N_;
  int p_;
  
  /*
   Enum for the units of energy
   */
  
  shared_ptr<VecOfMats<cmplx>::type> _gamma_;
  
  /*
   Calculate H vector (eq 42 and 43 in Lotan 2006)
   */
  VecOfMats<cmplx>::type calc_H(int i);
  
public:
  
  TorqueCalc(SHCalc shCalc, MyVector<VecOfMats<cmplx>::type> gradL,
             Constants consts, System sys, VecOfMats<cmplx>::type gamma,
             int p);
  
  void calc_tau();  // fill tau_
  
  /*
   Calculate inner product of two matrices as defined in equation 29 of Lotan
   2006
   */
  cmplx lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V, int p)
  {
    cmplx ip;
    int n, m;
    for (n = 0; n < p; n++)
    {
      for (m = -n; m <= -n; m++)
      {
        ip += U(n, m+p) * conj(V(n, m+p));
      }
    }
    return ip;
  }
  
};

#endif /* EnergyForce_h */

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
  int N_;  // number of molecules
  int p_;  // max number of poles
  shared_ptr<Constants> _const_;

  // result of energy calculation, internal units
  shared_ptr<MyVector<double> > _omega_;
public:
  EnergyCalc() { }
  
  EnergyCalc(shared_ptr<VecOfMats<cmplx>::type> _A,
             shared_ptr<VecOfMats<cmplx>::type> _L,
             shared_ptr<Constants> _const, int N, int p);
  
  EnergyCalc(shared_ptr<ASolver> _asolv);
  
//  EnergyCalc(ASolver asolv, Constants consts, int p);
  
  // fill omega_
  void calc_energy();
  
  // calculate energy of one molecule
  double calc_ei(int i);
  
  // get the energy for a specific molecule:
  double get_omega_i_int(int i)  { return _omega_->operator[](i); }
  // get all energy:
  shared_ptr<MyVector<double> > get_omega_int() { return _omega_; }
  
  // energy in kCal/mol:
  double get_omega_i_kcal(int i)
  { return _const_->convert_int_to_kcal_mol(_omega_->operator[](i)); }
  
  MyVector<double> get_omega_kcal()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = _const_->convert_int_to_kcal_mol(_omega_->operator[](n));
    return omeg;
  }
  
  // energy in kT:
  double get_omega_i_kT(int i)
  { return _const_->convert_int_to_kT(_omega_->operator[](i)); }
  MyVector<double> get_omega_kT()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = _const_->convert_int_to_kT(_omega_->operator[](n));
    return omeg;
  }
  
  // energy in joules/mol:
  double get_omega_i_jmol(int i)
  { return _const_->convert_int_to_jmol(_omega_->operator[](i)); }
  MyVector<double> get_omega_jmol()
  {
    MyVector<double> omeg(N_);
    for (int n = 0; n < N_; n++)
      omeg[n] = _const_->convert_int_to_jmol(_omega_->operator[](n));
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
  shared_ptr<Constants> _const_;
  
  shared_ptr<VecOfVecs<double>::type> _F_;
  
public:
  ForceCalc() { }
  
  ForceCalc(shared_ptr<VecOfMats<cmplx>::type> _A,
            shared_ptr<MyMatrix<VecOfMats<cmplx>::type > > _gradA,
            shared_ptr<VecOfMats<cmplx>::type> _L,
            shared_ptr<MyVector<VecOfMats<cmplx>::type > > _gradL,
            shared_ptr<Constants> con, int N, int p);
  
  ForceCalc(shared_ptr<ASolver> _asolv);
  
  void calc_force();  // fill F_
  
  // calculate force on one molecule
  MyVector<double> calc_fi(int i);
  
  MyVector<double> get_fi(int i)     { return _F_->operator[](i); }
  shared_ptr<VecOfVecs<double>::type> get_F()    { return _F_; }
  
};

/*
 Class for calculating the torque on every molecule in the system
 */
class TorqueCalc
{
protected:
  
  // outer vector has an entry for every molecule. Inner vector is the torque
  // on that molecule
  shared_ptr<VecOfVecs<double>::type> _tau_;
  
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr<BesselCalc> _bCalc_;
  shared_ptr< MyVector<VecOfMats<cmplx>::type > > _gradL_;
  
  shared_ptr<Constants> _consts_;
  shared_ptr<System> _sys_;
  shared_ptr<VecOfMats<cmplx>::type> _gamma_;
  
//  double epsS_;
  int N_;
  int p_;

//  /*
//   Calculate H vector (eq 42 and 43 in Lotan 2006)
//   */
//  VecOfMats<cmplx>::type calc_H(int i);
  
public:
  TorqueCalc() { }
  
  TorqueCalc(shared_ptr<SHCalc> _shCalc,
             shared_ptr<BesselCalc> _bCalc,
             shared_ptr<MyVector<VecOfMats<cmplx>::type> > _gradL,
             shared_ptr<VecOfMats<cmplx>::type> _gamma,
             shared_ptr<Constants> _consts,
             shared_ptr<System> sys, int p);
  
  TorqueCalc(shared_ptr<ASolver> _asolv);
  
  void calc_tau();  // fill tau_
  
  // calculate torque on one molecule
  MyVector<double> calc_tau_i(int i);
  
  /*
   Calculate H vector (eq 42 and 43 in Lotan 2006)
   */
  VecOfMats<cmplx>::type calc_H(int i);
  
  MyVector<double> get_taui(int i)     { return _tau_->operator[](i); }
  shared_ptr<VecOfVecs<double>::type > get_Tau()    { return _tau_; }
  
  /*
   Calculate inner product of two matrices as defined in equation 29 of Lotan
   2006
   */
  double lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V, int p)
  {
    double ip = 0;
    int n, m, mT;
    for (n = 0; n < p; n++)
    {
      for (m = -n; m <= n; m++)
      {
        mT = (m < 0) ? -1*m : m;
        ip += U(n, mT+p_).real()*V(n, mT+p_).real()
               + U(n, mT+p_).imag()*V(n, mT+p_).imag();
      }
    }
    return ip;
  }
  
};

/*
 Class for calculating energy force and torque in one place
 */
class PhysCalc
{
protected:
  int N_; // number of particles
  double unit_conv_; // Conversion factor for units
  string unit_; // String of the type of units
  vector<Pt> mol_pos_;
  
  shared_ptr<EnergyCalc> _eCalc_;
  shared_ptr<ForceCalc> _fCalc_;
  shared_ptr<TorqueCalc> _torCalc_;
  
  void compute_units( shared_ptr<Constants> cst, Units unit);
  
public:
  
  // constructor just requires an asolver
  PhysCalc(shared_ptr<ASolver> _asolv, Units unit = INTERNAL);
  
  MyVector<double> calc_force_i(int i)  { return _fCalc_->calc_fi(i); }
  MyVector<double> calc_tau_i(int i)    { return _torCalc_->calc_tau_i(i); }
  MyVector<double> calc_ei(int i)       { return _eCalc_->calc_ei(i); }
  
  void calc_force()   { _fCalc_->calc_force(); }
  void calc_energy()  { _eCalc_->calc_energy(); }
  void calc_torque()  { _torCalc_->calc_tau(); }
  void calc_all()     { calc_energy(); calc_force(); calc_torque(); }
  
  void print_all();
  
  shared_ptr<VecOfVecs<double>::type > get_Tau() { return _torCalc_->get_Tau();}
  shared_ptr<VecOfVecs<double>::type> get_F()    { return _fCalc_->get_F();}
  shared_ptr<MyVector<double> > get_omega() {return _eCalc_->get_omega_int();}
  
  MyVector<double> get_taui(int i) { return _torCalc_->get_taui(i); }
  MyVector<double> get_forcei(int i) { return _fCalc_->get_fi(i); }
  double get_omegai(int i) {return _eCalc_->get_omega_i_int(i);}
  
  MyVector<double> get_taui_conv(int i)
  { return _torCalc_->get_taui(i)*unit_conv_; }
  MyVector<double> get_forcei_conv(int i)
  { return _fCalc_->get_fi(i)*unit_conv_; }
  double get_omegai_conv(int i)
  {return _eCalc_->get_omega_i_int(i)*unit_conv_;}

};

#endif /* EnergyForce_h */

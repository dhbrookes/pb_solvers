//
//  BasePhysCalc.h
//  pbsam_xcode
//
//  Created by David Brookes on 9/1/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BasePhysCalc_h
#define BasePhysCalc_h

#include <stdio.h>
#include <memory>
#include <vector>
#include "util.h"
#include "BaseSys.h"

using namespace std;
/*
 Base class for energy calculations
 */
class BaseEnergyCalc
{
protected:
  shared_ptr<vector<double> > omega_;
  
public:
  BaseEnergyCalc() { }
  BaseEnergyCalc(int N) { omega_ = make_shared<vector<double> > (N); }
  
  shared_ptr<vector<double> > get_omega() const { return omega_; }
  const double get_ei(int i) const              { return (*omega_)[i]; }
  
  virtual void calc_all_energy() { }
};

/*
 Base class for force calculations
 */
class BaseForceCalc
{
protected:
  shared_ptr<vector<Pt> > _F_;
  
public:
  BaseForceCalc() { }
  BaseForceCalc(int N) { _F_ = make_shared<vector<Pt> >(N); }
  
  shared_ptr<vector<Pt> > get_F() const { return _F_; }
  Pt get_forcei(int i) const                    { return (*_F_)[i]; }
  
  virtual void calc_all_f() { }
};


class BaseTorqueCalc
{
protected:
  shared_ptr<vector<Pt> > _tau_;
  
public:
  BaseTorqueCalc() { }
  BaseTorqueCalc(int N) { _tau_ = make_shared<vector<Pt> > (N); }
  
  shared_ptr<vector<Pt> > get_tau() const   { return _tau_; }
  Pt get_taui(int i) const                  { return (*_tau_)[i]; }
  
  virtual void calc_all_tau() { }
};

/*
 Base class for calculations of physical quantities
 */
class BasePhysCalc
{
protected:
  int N_; // number of particles
  double unit_conv_; // Conversion factor for units
  string unit_; // String of the type of units
  string outfname_; // where you want the info printed to
  
  void compute_units( shared_ptr<Constants> cst, Units unit)
  {
    if (unit==INTERNAL)
    {
      unit_ = "Internal";
      unit_conv_ = 1.0;
    } else if (unit == KCALMOL)
    {
      unit_  = "kCal/Mol";
      unit_conv_ = cst->convert_int_to_kcal_mol(1.0);
    } else if (unit == JMOL)
    {
      unit_  = "Joules/Mol";
      unit_conv_ = cst->convert_int_to_jmol(1.0);
    } else if (unit == kT)
    {
      unit_  = "kT";
      unit_conv_ = cst->convert_int_to_kT(1.0);
    }
  }
  
public:
  BasePhysCalc() { }
  
  BasePhysCalc(int N, shared_ptr<Constants> cst,
               string outfname, Units units=INTERNAL)
  : N_(N), outfname_(outfname)
  {
    compute_units(cst, units);
  }
  
  virtual void calc_force() {}
  virtual void calc_energy() {}
  virtual void calc_torque() {}
  
  
  void calc_all()     { calc_energy(); calc_force(); calc_torque(); }
  
  virtual void print_all() { }
  
  virtual shared_ptr<vector<Pt> > get_Tau()
  { return make_shared<vector<Pt> > ();  }
  virtual shared_ptr<vector<Pt> > get_F()
  { return make_shared<vector<Pt> > ();   }
  virtual shared_ptr<vector<double> > get_omega()
  { return make_shared<vector<double> > (); }
  
  virtual Pt get_taui(int i) { return Pt(); }
  virtual Pt get_forcei(int i) { return Pt (); }
  virtual double get_omegai(int i) {return 0; }
  virtual double calc_ei(int i) {return 0; }
  
  virtual Pt get_moli_pos(int i) { return Pt(); }
};

#endif /* BasePhysCalc_h */

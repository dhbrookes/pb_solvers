//
//  PhysCalcCalc.hpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef PhysCalc_h
#define PhysCalc_h

#include <stdio.h>
#include <memory>
#include "Solver.h"
#include <map>


/*
 Class for calculating interaction energy of a MoleculeSAM given H^(I,k)
 and LHN^(I,k) matrices
 */
class EnergyCalc
{
  shared_ptr<vector<double> > omega_;
  
public:
  EnergyCalc(int I) { omega_ = make_shared<vector<double> > (I); }
  
  double calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN);
  void calc_all_energy(vector<shared_ptr<HMatrix> > H,
                       vector<shared_ptr<LHNMatrix> > LHN);
  
  shared_ptr<vector<double> > get_omega() { return omega_; }
};



class ForceCalc
{
protected:
  shared_ptr<SHCalc> shcalc_;
  shared_ptr<BesselCalc> bcalc_;
  
  vector<vector<Pt> > forces_;
  shared_ptr<vector<Pt> > totForces_; // forces on each molecule
  int I_; // number of molecules in the system
  vector<int> ks_;
  
  double eps_s_;
  
public:
  
  // number of molecules, Ik of all molecules, solvent dielectric:
  ForceCalc(int I, vector<int> ki, double es, shared_ptr<SHCalc> shcalc,
            shared_ptr<BesselCalc> bcalc)
  :shcalc_(shcalc), bcalc_(bcalc), forces_(I), I_(I), ks_(ki), eps_s_(es)
  {
    totForces_ = make_shared<vector<Pt> >(I_);
  }
  
  /*
   Calculate translational force on a MoleculeSAM given dI_LHN^(I,k), H^(I,k),
   LHN^(I,k) and dI_H^(I,k)
   */
  void calc_fI(int I, shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
               shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
  
  void calc_all_f(vector<shared_ptr<HMatrix> > H,
                  vector<shared_ptr<LHNMatrix> > LHN,
                  vector<vector<shared_ptr<GradHMatrix> > > dH,
                  vector<vector<shared_ptr<GradLHNMatrix> > > dLHN);
  
  //calc force at a point
//  Ptx calc_fp(Pt P, shared_ptr<BaseMolecule> mol,
//              shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
//              shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
  
  shared_ptr<vector<Pt> > get_all_f()
  {
    return totForces_;
  }
  
  vector<Pt> get_all_fIk(int I)  {return forces_[I];}
};


class TorqueCalc
{
protected:
  shared_ptr<vector<Pt> > torques_;
  int I_; // number of MoleculeSAMs in the system
  
public:
  TorqueCalc(int I) : I_(I)  { torques_ = make_shared<vector<Pt> >(I_); }
  
  void calc_all_tau(shared_ptr<System> sys, shared_ptr<ForceCalc> fcalc);
  Pt calc_tauI(int i, shared_ptr<BaseMolecule> mol, shared_ptr<ForceCalc> fcalc);
  
  Pt cross_prod(Pt a, Pt b);
  
  shared_ptr<vector<Pt> > get_all_tau()   { return torques_;}
};


/*
 Base class for calculations of physical quantities
 */
class BasePhysCalc
{
public:
  BasePhysCalc() { }
  
  virtual void calc_force() {}
  virtual void calc_energy() {}
  virtual void calc_torque() {}
  
  virtual void print_all() { }
  
  virtual shared_ptr<vector<Pt> > get_Tau()
  { return make_shared<vector<Pt> > ();  }
  virtual shared_ptr<vector<Pt> > get_F()
  { return make_shared<vector<Pt> > ();   }
  virtual shared_ptr<vector<double> > get_omega()
  { return make_shared<vector<double> > (); }
  
};

/*
 Class for calculating energy force and torque in one place
 */
class PhysCalc : public BasePhysCalc
{
protected:
  int N_; // number of particles
  double unit_conv_; // Conversion factor for units
  string unit_; // String of the type of units
  shared_ptr<System> _sys_; // System
  
  shared_ptr<EnergyCalc> _eCalc_;
  shared_ptr<ForceCalc> _fCalc_;
  shared_ptr<TorqueCalc> _torCalc_;
  
  shared_ptr<Solver> _solv_;
  shared_ptr<GradSolver> _gradSolv_;
  
  string outfname_; // where you want the info printed to
  
  void compute_units( shared_ptr<Constants> cst, Units unit);
  
public:
  
  // constructor just requires an asolver
  PhysCalc(shared_ptr<Solver> _solv, shared_ptr<GradSolver> _gradsolv,
           string outfname, Units unit = INTERNAL);
  
//  void calc_force_interact()   { _fCalc_->calc_force_interact(_sys_); }
  void calc_force();
  
  
  void calc_energy();
  void calc_torque();
  void calc_all()     { calc_energy(); calc_force(); calc_torque(); }
  
//  void print_all();
  
  shared_ptr<vector<Pt> > get_Tau() { return _torCalc_->get_all_tau(); }
  shared_ptr<vector<Pt> > get_F() { return _fCalc_->get_all_f(); }
  shared_ptr<vector<double> > get_omega() { return _eCalc_->get_omega(); }
  
  Pt get_moli_pos( int i) { return _sys_->get_cogi(i); }
  
};

#endif /* PhysCalcCalc_h */





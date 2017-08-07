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
#include "BasePhysCalc.h"
#include <map>

/*
 Class for calculating interaction energy of a MoleculeSAM given H^(I,k)
 and LHN^(I,k) matrices
 */
class EnergyCalcSAM : public BaseEnergyCalc
{
  double epsS;

public:
  EnergyCalcSAM(int I, double eps_s) : BaseEnergyCalc(I), epsS(eps_s) { }

  double calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN);

  void calc_all_energy(vector<shared_ptr<HMatrix> > H,
                       vector<shared_ptr<LHNMatrix> > LHN);

};

class ForceCalcSAM : public BaseForceCalc
{
protected:
  shared_ptr<SHCalc> shcalc_;
  shared_ptr<BesselCalc> bcalc_;

  // Force on each CG sphere of each molecule
  vector<vector<Pt> > forces_;

  int I_; // number of molecules in the system
  vector<int> ks_;

  double eps_s_;

  /*
   Calculate translational force on a MoleculeSAM given dI_LHN^(I,k), H^(I,k),
   LHN^(I,k) and dI_H^(I,k)
   */
  void calc_fI(int I, shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
               shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);

public:

  // number of molecules, Ik of all molecules, solvent dielectric:
  ForceCalcSAM(int I, vector<int> ki, double es, shared_ptr<SHCalc> shcalc,
            shared_ptr<BesselCalc> bcalc)
  :BaseForceCalc(I), shcalc_(shcalc), bcalc_(bcalc), I_(I), ks_(ki), eps_s_(es),
  forces_(I)
  {
  }

  void calc_all_f(vector<shared_ptr<HMatrix> > H,
                  vector<shared_ptr<LHNMatrix> > LHN,
                  vector<vector<shared_ptr<GradHMatrix> > > dH,
                  vector<vector<shared_ptr<GradLHNMatrix> > > dLHN);

  vector<Pt> get_all_fIk(int I)  {return forces_[I];}
};


class TorqueCalcSAM : public BaseTorqueCalc
{
protected:
  int I_;

public:
  TorqueCalcSAM(int I) : BaseTorqueCalc(I), I_(I)  { }

  void calc_all_tau(shared_ptr<SystemSAM> sys, shared_ptr<ForceCalcSAM> fcalc);
  Pt calc_tauI(int i, shared_ptr<BaseMolecule> mol, shared_ptr<ForceCalcSAM> fcalc);

  Pt cross_prod(Pt a, Pt b);
};

/*
 Class for calculating energy force and torque in one place
 */
class PhysCalcSAM : public BasePhysCalc
{
protected:
  shared_ptr<Solver> _solv_;
  shared_ptr<GradSolver> _gradSolv_;

  shared_ptr<EnergyCalcSAM> _eCalc_;
  shared_ptr<ForceCalcSAM>  _fCalc_;
  shared_ptr<TorqueCalcSAM> _torCalc_;

  shared_ptr<SystemSAM> _sys_;

public:

  // constructor just requires an asolver
  PhysCalcSAM(shared_ptr<Solver> _solv, shared_ptr<GradSolver> _gradsolv,
           string outfname, Units unit = INTERNAL);

  void calc_force();
  void calc_energy();
  void calc_torque();

  void print_all();

  shared_ptr<vector<Pt> > get_Tau() { return _torCalc_->get_tau(); }
  shared_ptr<vector<Pt> > get_F()   { return _fCalc_->get_F(); }
  shared_ptr<vector<double> > get_omega() { return _eCalc_->get_omega(); }

  Pt get_taui(int i)        { return _torCalc_->get_taui(i); }
  Pt get_forcei(int i)      { return _fCalc_->get_forcei(i); }
  double get_omegai(int i)  { return _eCalc_->get_ei(i); }
  Pt get_taui_conv(int i)        { return _torCalc_->get_taui(i)*unit_conv_; }
  Pt get_forcei_conv(int i)      { return _fCalc_->get_forcei(i)*unit_conv_; }
  double get_omegai_conv(int i)  { return _eCalc_->get_ei(i)*unit_conv_; }
  Pt get_moli_pos(int i)    { return _sys_->get_cogi(i); }

};

#endif /* PhysCalcCalc_h */




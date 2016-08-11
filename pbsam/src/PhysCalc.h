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


/*
 Class for calculating interaction energy of a MoleculeSAM given H^(I,k)
 and LHN^(I,k) matrices
 */
class EnergyCalc
{
public:
  EnergyCalc() { }
  
  double calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN);
  vector<double> calc_all_energy(vector<shared_ptr<HMatrix> > H,
                                 vector<shared_ptr<LHNMatrix> > LHN);
};



class ForceCalc
{
protected:
  shared_ptr<SHCalc> shcalc_;
  shared_ptr<BesselCalc> bcalc_;
  
  vector<vector<Pt> > forces_;
  int I_; // number of MoleculeSAMs in the system
  vector<int> ks_;
  
  double eps_s_;
  
public:
  ForceCalc(int I, vector<int> ki, double es, shared_ptr<SHCalc> shcalc,
            shared_ptr<BesselCalc> bcalc)
  :shcalc_(shcalc), bcalc_(bcalc), forces_(I), I_(I), ks_(ki), eps_s_(es)  { }
  
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
  
  vector<Pt> get_all_f()
  {
    vector<Pt> forc(I_, 0.0);
    for (int i=0; i<I_; i++)
      for (int k=0; k<ks_[i]; k++)
        forc[i] += forces_[i][k];
    
    return forc;
  }
  
  vector<Pt> get_all_fIk(int I)  {return forces_[I];}
};


class TorqueCalc
{
protected:
  vector<Pt> torques_;
  int I_; // number of MoleculeSAMs in the system
  
public:
  TorqueCalc(int I) :  torques_(I), I_(I)  { }
  
  void calc_all_tau(shared_ptr<System> sys, shared_ptr<ForceCalc> fcalc);
  Pt calc_tauI(int i, shared_ptr<BaseMolecule> mol, shared_ptr<ForceCalc> fcalc);
  
  Pt cross_prod(Pt a, Pt b);
  
  vector<Pt> get_all_tau()   {return torques_;}
};

#endif /* PhysCalcCalc_h */

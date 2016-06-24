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
#include "Solvmat.h"
#include "Gradsolvmat.h"


/*
 Class for calculating interaction energy of a molecule given H^(I,k)
 and LHN^(I,k) matrices
 */
class EnergyCalc
{
public:
  EnergyCalc() { }
  
  cmplx calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN);
  
};



class ForceCalc
{
protected:
  shared_ptr<SHCalc> shcalc_;
  shared_ptr<BesselCalc> bcalc_;
  
public:
  ForceCalc(shared_ptr<SHCalc> shcalc, shared_ptr<BesselCalc> bcalc)
  :shcalc_(shcalc), bcalc_(bcalc)
  {
  }
  
  /*
   Calculate translational force on a molecule given dI_LHN^(I,k), H^(I,k),
   LHN^(I,k) and dI_H^(I,k)
   */
  Ptx calc_fI(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
              shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
  
  //calc force at a point
  Ptx calc_fp(Pt P, shared_ptr<Molecule> mol,
              shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
              shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
};


class TorqueCalc
{
public:
  TorqueCalc()  { }
  
  Ptx calc_tauI(shared_ptr<Molecule> mol, shared_ptr<ForceCalc> fcalc,
                shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
                shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
  
  Ptx cross_prod(Pt a, Ptx b);
};

#endif /* PhysCalcCalc_h */

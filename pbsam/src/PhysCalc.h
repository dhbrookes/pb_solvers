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


/*
 Calculate translational force on a molecule given dI_LHN^(I,k), H^(I,k),
 LHN^(I,k) and dI_H^(I,k)
 */
class ForceCalc
{
public:
  ForceCalc() { }
  
  Ptx calc_force(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
                 shared_ptr<GradHMatrix> dH, shared_ptr<GradLHNMatrix> dLHN);
};

#endif /* PhysCalcCalc_h */

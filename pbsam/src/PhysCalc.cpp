//
//  PhysCalc.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "PhysCalc.h"


cmplx EnergyCalc::calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN)
{
  cmplx E = 0;
  for (int k = 0; k < H->get_ns(); k++)
    for (int n = 0; n < H->get_p(); n++)
      for (int m = - n; m < n+1; m++)
        E += LHN->get_mat_knm(k, n, m) * conj(H->get_mat_knm(k, n, m));
  
  return E;
}


Ptx ForceCalc::calc_force(shared_ptr<HMatrix> H,
                          shared_ptr<LHNMatrix> LHN,
                          shared_ptr<GradHMatrix> dH,
                          shared_ptr<GradLHNMatrix> dLHN)
{
  Ptx f, fIk, inner1, inner2;
  for (int k = 0; k < H->get_ns(); k++)
    for (int n = 0; n < H->get_p(); n++)
      for (int m = - n; m < n+1; m++)
      {
        inner1 = dLHN->get_mat_knm(k, n, m) * H->get_mat_knm(k, n, m);
        inner2 = H->get_mat_knm(k, n, m) * LHN->get_mat_knm(k, n, m);
        fIk = inner1 * (-1) + inner2;
        fIk *= -1;
        f += fIk;
      }
  return f;
}
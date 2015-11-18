//
//  EnergForce.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 11/18/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "EnergyForce.h"


EnergyCalc::EnergyCalc(VecOfMats<cmplx>::type A, VecOfMats<cmplx>::type L,
                       double epsS, int N, int p)
:N_(N), epsS_(epsS), omega_(N), p_(p)
{
  _A_ = make_shared<VecOfMats<cmplx>::type> (A);
  _L_ = make_shared<VecOfMats<cmplx>::type> (L);
  
  calc_energy();
}

void EnergyCalc::calc_energy()
{
  int i, n, m;
  cmplx ip, unm, vnm;
  MyMatrix<cmplx> Li, Ai;
  for (i = 0; i < N_; i++)
  {
    Li = _L_->operator[](i);
    Ai = _A_->operator[](i);
    
    // calculate inner product (as defined in eq 29 of Lotan 2006):
    ip = cmplx (0, 0);
    for (n = 0; n < p_; n++)
    {
      for (m = -n; m < p_; m++)
      {
        unm = Li(n, m + p_);
        vnm = Ai(n, m + p_);
        ip += unm * conj(vnm);
      }
    }
    omega_.set_val(i, ip * (1/epsS_));
  }
}


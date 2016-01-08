//
//  EnergForce.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 11/18/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "EnergyForce.h"

EnergyCalc::EnergyCalc(VecOfMats<cmplx>::type A, VecOfMats<cmplx>::type L,
                       Constants con, int N, int p)
:N_(N), const_(con), omega_(N), p_(p)
{
  _A_ = make_shared<VecOfMats<cmplx>::type> (A);
  _L_ = make_shared<VecOfMats<cmplx>::type> (L);
  
  calc_energy();
}

void EnergyCalc::calc_energy()
{
  int i, n, m; cmplx unm, vnm; double ip; MyMatrix<cmplx> Li, Ai;
  
  for (i = 0; i < N_; i++)
  {
    Li = _L_->operator[](i);
    Ai = _A_->operator[](i);
    
    // calculate inner product (as defined in eq 29 of Lotan 2006):
    ip = 0.0;
    for (n = 0; n < p_; n++)
    {
      for (m = -n; m <= n; m++)
      {
        unm = Li(n, m + p_);
        vnm = Ai(n, m + p_);
        ip += unm.real()*vnm.real() + unm.imag()*vnm.imag();
      }
    }
    omega_.set_val(i, ip * (1/const_.get_dielectric_water()));
  }
}


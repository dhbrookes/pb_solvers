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


ForceCalc::ForceCalc(VecOfMats<cmplx>::type A,
                     MyMatrix<VecOfMats<cmplx>::type > gradA_,
                     VecOfMats<cmplx>::type L,
                     MyVector<VecOfMats<cmplx>::type > gradL_,
                     double epsS, int N, int p)
:N_(N), epsS_(epsS), F_(N), p_(p)
{
  _A_ = make_shared<VecOfMats<cmplx>::type> (A);
  _L_ = make_shared<VecOfMats<cmplx>::type> (L);
  
  _gradA_ = make_shared<MyMatrix<VecOfMats<cmplx>::type> > (gradA);
  _gradL_ = make_shared<MyVector<VecOfMats<cmplx>::type> > (gradL);
  
}


void ForceCalc::calc_force()
{
  int i, j, n, m;
  cmplx ip1, ip2, unm1, vnm1, unm2, vnm2;
  VecOfMats<cmplx> gLi, gAi;
  MyMatrix<cmplx> Li, Ai;
  MyVector<cmplx> inner;
  for (i = 0; i < N_; i++)
  {
    Li = _L_->operator[](i);
    Ai = _A_->operator[](i);
    
    gLi = _gradL_->operator[](i)
    gAi = _gradA_->operator[](i);
    inner = MyVector<cmplx> (3);
    for (j = 0; j < 3; j++)  // for each component of the gradient
    {
      ip1 = cmplx(0, 0);
      ip2 = cmplx(0, 0);
      for (n = 0; n < p_; n++)
      {
        for (m = -n; m < p_; m++)
        {
          unm1 = gLi[j](n, m + p_);
          vnm1 = Ai(n, m + p_);
          ip1 += unm1 * conj(vnm1);
          
          unm2 = Li(n, m + p);
          vnm2 = gAi[j](n, m + p);
          ip2 += unm1 * conj(vnm1);
        }
      }
      inner.set_val(j, 1/epsS_ * (ip1 + ip2))
    }
    F_.set_val(i, inner);
  }
}
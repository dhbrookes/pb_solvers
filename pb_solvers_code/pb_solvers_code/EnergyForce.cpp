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

ForceCalc::ForceCalc(VecOfMats<cmplx>::type A,
                     MyMatrix<VecOfMats<cmplx>::type > gradA_,
                     VecOfMats<cmplx>::type L,
                     MyVector<VecOfMats<cmplx>::type > gradL_,
                     Constants con, int N, int p)
:N_(N), const_(con), F_(N), p_(p)
{
  _A_ = make_shared<VecOfMats<cmplx>::type> (A);
  _L_ = make_shared<VecOfMats<cmplx>::type> (L);
  
  _gradA_ = make_shared<MyMatrix<VecOfMats<cmplx>::type> > (gradA_);
  _gradL_ = make_shared<MyVector<VecOfMats<cmplx>::type> > (gradL_);
  
  calc_force();
}


void ForceCalc::calc_force()
{
  int i, j, n, m;
  cmplx unm1, vnm1, unm2, vnm2;
  double ip1, ip2;
  VecOfMats<cmplx>::type gLi, gAi;
  MyMatrix<cmplx> Li, Ai;
  MyVector<double> inner;
  for (i = 0; i < N_; i++)
  {
    Li = _L_->operator[](i);
    Ai = _A_->operator[](i);
    
    gLi = _gradL_->operator[](i);
    gAi = _gradA_->operator()(i, i);
    inner = MyVector<double> (3);
    for (j = 0; j < 3; j++)  // for each component of the gradient
    {
      ip1 = 0.0;
      ip2 = 0.0;
      for (n = 0; n < p_; n++)
      {
        for (m = -n; m < p_; m++)
        {
          unm1 = gLi[j](n, m + p_);
          vnm1 = Ai(n, m + p_);
          ip1 += unm1.real()*vnm1.real() + unm1.imag()*vnm1.imag();
          
          unm2 = Li(n, m + p_);
          vnm2 = gAi[j](n, m + p_);
          ip2 += unm2.real()*vnm2.real() + unm2.imag()*vnm2.imag();
        }
      }
      inner.set_val(j, -1.0/const_.get_dielectric_water() * (ip1 + ip2));
    }
    F_.set_val(i, inner);
  }
}
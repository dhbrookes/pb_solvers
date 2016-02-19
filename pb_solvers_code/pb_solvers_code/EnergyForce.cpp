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
        for (m = -n; m <= n; m++)
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

TorqueCalc::TorqueCalc(SHCalc shCalc, BesselCalc bCalc,
                       MyVector<VecOfMats<cmplx>::type> gradL,
                       VecOfMats<cmplx>::type gamma,
                       Constants consts, System sys, int p)
: N_(sys.get_n()), p_(p), tau_(sys.get_n()), consts_(consts),
epsS_(consts.get_dielectric_water())
{
  _shCalc_ = make_shared<SHCalc> (shCalc);
  _bCalc_  = make_shared<BesselCalc> (bCalc);
  _gradL_  = make_shared<MyVector<VecOfMats<cmplx>::type > > (gradL);
  _gamma_ = make_shared<VecOfMats<cmplx>::type > (gamma);
  _sys_    = make_shared<System> (sys);
  
  calc_tau();
}


VecOfMats<cmplx>::type TorqueCalc::calc_H(int i)
{
  VecOfMats<cmplx>::type H (3);
  int mi = _sys_->get_Mi(i);
  cmplx sh, h, gam;
  double scale, qij;
  double lambda  = _sys_->get_lambda();
  MyMatrix<cmplx> gamma_i = _gamma_->operator[](i);
  vector<double> bessI;
  MyMatrix<cmplx> Hx (p_, 2*p_), Hy (p_, 2*p_), Hz (p_, 2*p_);
  
  int j, n, m;
  for (j = 0; j < mi; j++)
  {
    qij = _sys_->get_qij(i, j);
    Pt pt = _sys_->get_posij(i, j);
    _shCalc_->calc_sh(pt.theta(),pt.phi());
    scale = 1.0;
    
    if (_sys_->get_ai(i) == 0)
      bessI = _bCalc_->calc_mbfI(p_, _sys_->get_consts().get_kappa()*pt.r());
    else
      bessI = _bCalc_->calc_mbfI(p_, 0.0);

    for (n = 0; n < p_; n++)
    {
      gam = gamma_i(n, n);
      for (m = 0; m <= n; m++)
      {
        sh = _shCalc_->get_result(n, m);
        h = bessI[n] * qij * scale * sh * gam;
        Hx(n, m+p_) += h * pt.x();
        Hy(n, m+p_) += h * pt.y();
        Hz(n, m+p_) += h * pt.z();
      }
      scale *= (pt.r()/lambda);
    }
  }
  H.set_val(0, Hx);
  H.set_val(1, Hy);
  H.set_val(2, Hz);
  
  return H;
}


void TorqueCalc::calc_tau()
{
  MyVector<double> tau_i;
  VecOfMats<cmplx>::type Hi;
  VecOfMats<cmplx>::type gLi;
  for (int i = 0; i < N_; i++)
  {
    Hi    = calc_H(i);
    tau_i = MyVector<double> (3);
    gLi   = _gradL_-> operator[](i);
    
    //perform cross product:
    tau_i.set_val(0, 1/epsS_ * (lotan_inner_prod(gLi[1], Hi[2], p_)
                  - lotan_inner_prod(gLi[2], Hi[1], p_)));
    
    tau_i.set_val(1, 1/epsS_ * (lotan_inner_prod(gLi[2], Hi[0], p_)
                  - lotan_inner_prod(gLi[0], Hi[2], p_)));
    
    tau_i.set_val(2, 1/epsS_ * (lotan_inner_prod(gLi[0], Hi[1], p_)
                  - lotan_inner_prod(gLi[1], Hi[0], p_)));
    
    tau_.set_val(i, tau_i);
  }
}


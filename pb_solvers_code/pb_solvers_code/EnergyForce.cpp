//
//  EnergForce.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 11/18/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "EnergyForce.h"

EnergyCalc::EnergyCalc(shared_ptr<VecOfMats<cmplx>::type> _A,
                       shared_ptr<VecOfMats<cmplx>::type> _L,
                       shared_ptr<Constants> _const, int N, int p)
:N_(N), _const_(_const), p_(p), _A_(_A), _L_(_L)
{
  _omega_ = make_shared<MyVector<double> > (N_);
  calc_energy();
}

EnergyCalc::EnergyCalc(shared_ptr<ASolver> _asolv)
:_A_(_asolv->get_A()), _L_(_asolv->get_L()), _const_(_asolv->get_consts()),
N_(_asolv->get_N()), p_(_asolv->get_p())
{
  _omega_ = make_shared<MyVector<double> > (N_);
//  calc_energy();
}

double EnergyCalc::calc_ei(int i)
{
  double ei;
  int n, m;
  cmplx unm, vnm;
  MyMatrix<cmplx> Li = _L_->operator[](i);
  MyMatrix<cmplx> Ai = _A_->operator[](i);
  
  // calculate inner product (as defined in eq 29 of Lotan 2006):
  ei = 0.0;
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      unm = Li(n, m + p_);
      vnm = Ai(n, m + p_);
      ei += unm.real()*vnm.real() + unm.imag()*vnm.imag();
    }
  }
  ei *= (1/_const_->get_dielectric_water());
  return ei;
}

void EnergyCalc::calc_energy()
{
  for (int i = 0; i < N_; i++)
  {
    _omega_->set_val(i, calc_ei(i));
  }
}

ForceCalc::ForceCalc(shared_ptr<VecOfMats<cmplx>::type> _A,
                     shared_ptr<MyMatrix<VecOfMats<cmplx>::type > > _gradA,
                     shared_ptr<VecOfMats<cmplx>::type> _L,
                     shared_ptr<MyVector<VecOfMats<cmplx>::type > > _gradL,
                     shared_ptr<Constants> _con, int N, int p)
:N_(N), _const_(_con), p_(p), _gradA_(_gradA), _A_(_A), _L_(_L),
_gradL_(_gradL)
{
  _F_ = make_shared<VecOfVecs<double>::type > (N_, MyVector<double> (3));
  calc_force();
}

ForceCalc::ForceCalc(shared_ptr<ASolver> _asolv)
:_A_(_asolv->get_A()), _gradA_(_asolv->get_gradA()), _L_(_asolv->get_L()),
_gradL_(_asolv->get_gradL()), _const_(_asolv->get_consts()),
N_(_asolv->get_N()), p_(_asolv->get_p())
{
  _F_ = make_shared<VecOfVecs<double>::type > (N_, MyVector<double> (3));
//  calc_force();
}

MyVector<double> ForceCalc::calc_fi(int i)
{
  int j, n, m;
  cmplx unm1, vnm1, unm2, vnm2;
  double ip1, ip2;
  VecOfMats<cmplx>::type gLi, gAi;
  MyMatrix<cmplx> Li, Ai;
  MyVector<double> fi;

  Li = _L_->operator[](i);
  Ai = _A_->operator[](i);
  
  gLi = _gradL_->operator[](i);
  gAi = _gradA_->operator()(i, i);
  fi = MyVector<double> (3);
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
    
    fi.set_val(j, -1.0/_const_->get_dielectric_water() * (ip1 + ip2));
  }
  return fi;

}


void ForceCalc::calc_force()
{
  for (int i = 0; i < N_; i++)
  {
    _F_->set_val(i, calc_fi(i));
  }
}

TorqueCalc::TorqueCalc(shared_ptr<SHCalc> _shCalc,
                       shared_ptr<BesselCalc> _bCalc,
                       shared_ptr<MyVector<VecOfMats<cmplx>::type> > _gradL,
                       shared_ptr<VecOfMats<cmplx>::type> _gamma,
                       shared_ptr<Constants> _consts,
                       shared_ptr<System> _sys, int p)
: N_(_sys->get_n()), p_(p), _consts_(_consts),
_shCalc_(_shCalc), _bCalc_(_bCalc), _gradL_(_gradL), _gamma_(_gamma)
{
  _tau_ = make_shared<VecOfVecs<double>::type > (N_, MyVector<double> (3));
  calc_tau();
}

TorqueCalc::TorqueCalc(shared_ptr<ASolver> _asolv)
:N_(_asolv->get_N()), p_(_asolv->get_p()),
_consts_(_asolv->get_consts()), _shCalc_(_asolv->get_sh()),
_bCalc_(_asolv->get_bessel()),
_gamma_(_asolv->get_gamma()),
_sys_(_asolv->get_sys()),
_gradL_(_asolv->get_gradL())
{
  _tau_ = make_shared<VecOfVecs<double>::type > (N_, MyVector<double> (3));
//  calc_tau();
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
      bessI = _bCalc_->calc_mbfI(p_, _consts_->get_kappa()*pt.r());
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


MyVector<double> TorqueCalc::calc_tau_i(int i)
{
  MyVector<double> tau_i;
  VecOfMats<cmplx>::type Hi;
  VecOfMats<cmplx>::type gLi;

  Hi    = calc_H(i);
  tau_i = MyVector<double> (3);
  gLi   = _gradL_-> operator[](i);
  
  //perform cross product:
  tau_i.set_val(0, 1/_consts_->get_dielectric_water()
                * (lotan_inner_prod(gLi[1], Hi[2], p_)
                   - lotan_inner_prod(gLi[2], Hi[1], p_)));
  
  tau_i.set_val(1, 1/_consts_->get_dielectric_water() *
                (lotan_inner_prod(gLi[2], Hi[0], p_)
                 - lotan_inner_prod(gLi[0], Hi[2], p_)));
  
  tau_i.set_val(2, 1/_consts_->get_dielectric_water() *
                (lotan_inner_prod(gLi[0], Hi[1], p_)
                 - lotan_inner_prod(gLi[1], Hi[0], p_)));
  return tau_i;
}

void TorqueCalc::calc_tau()
{
  for (int i = 0; i < N_; i++)
  {
    _tau_->set_val(i, calc_tau_i(i));
  }
}


PhysCalc::PhysCalc(shared_ptr<ASolver> _asolv)
{
  _eCalc_ = make_shared<EnergyCalc>(_asolv);
  _fCalc_ = make_shared<ForceCalc>(_asolv);
  _torCalc_ = make_shared<TorqueCalc>(_asolv);
}


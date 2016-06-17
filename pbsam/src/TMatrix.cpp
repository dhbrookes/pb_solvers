//
//  TMatrix.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/15/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "TMatrix.h"

TMatrix::TMatrix(int p, shared_ptr<System> _sys,
                 shared_ptr<SHCalc> _shcalc,
                 shared_ptr<Constants> _consts,
                 shared_ptr<BesselCalc> _besselcalc,
                 shared_ptr<ReExpCoeffsConstants> _reexpconsts)
:p_(p), kappa_(_consts->get_kappa()), Nmol_(_sys->get_n()), _system_(_sys),
_besselCalc_(_besselcalc), _shCalc_(_shcalc)
{
  int total_spheres=0;
  Nsi_ = vector<int> (Nmol_);
  for (int I = 0; I < Nmol_; I++)
  {
    total_spheres += _sys->get_Ns_i(I);
    Nsi_[I] = _sys->get_Ns_i(I);
  }
  T_.reserve(total_spheres*total_spheres);
  
  update_vals(_sys, _shcalc, _besselcalc, _reexpconsts);
}


void TMatrix::update_vals(shared_ptr<System> _sys, shared_ptr<SHCalc> _shcalc,
                          shared_ptr<BesselCalc> _besselcalc,
                          shared_ptr<ReExpCoeffsConstants> _reexpconsts)
{
  T_.clear();
  idxMap_.clear();
  
  int idx = 0;
  Pt c_Ik, c_Jl, v;
  vector<int> idx_vec;
  for (int I = 0; I < _sys->get_n(); I++)
  {
    for (int J = 0; J < _sys->get_n(); J++)
    {
      for (int k = 0; k < _sys->get_Ns_i(I); k++)
      {
        c_Ik = _sys->get_centerik(I, k);
        for (int l = 0; l < _sys->get_Ns_i(J); l++)
        {
          if (I==J && k==l) continue;
          
          idx_vec = {I, k, J, l};
          c_Jl = _sys->get_centerik(J, l);
          if (I==J && c_Ik.dist(c_Jl) < 5.0)
          {
            idxMap_[idx_vec] = -1;
            continue;
          }
          vector<double> besselK = _besselcalc->calc_mbfK(2*p_, kappa_ * v.r());
          v = _sys->get_pbc_dist_vec_base(c_Ik, c_Jl);
          _shcalc->calc_sh(v.theta(), v.phi());
          auto re_exp = make_shared<ReExpCoeffs>(p_, v,
                                                 _shcalc->get_full_result(),
                                                 besselK, _reexpconsts,
                                                 _sys->get_lambda(), false);
          T_.push_back(re_exp);
          idxMap_[idx_vec] = idx;
          idx++;
        }
      }
    }
  }
}

// Perform local expansion from J, l onto I, k
MyMatrix<cmplx> TMatrix::re_expandX_local(vector<vector<double> > X,
                                          int I, int k,
                                          int J, int l)
{
  int h, n, m;
  cmplx val;
  double chgscl, rscl;
  MyMatrix<cmplx> Z(p_, 2*p_+1);
  vector<int> exp_pts = _system_->get_gdpt_expij(J, l);
  for (h = 0; h < X[l].size(); h++)
  {
    Pt sph_dist = _system_->get_centerik(I, k) - _system_->get_centerik(J, l);
    Pt loc = _system_->get_gridijh(J, l, exp_pts[h]) - sph_dist;
//    cout << "Hval " << X[l][h] << " at: " << loc.x() << ", " << loc.y() << ", " << loc.z() << " rho: " << loc.r() << endl;
    _shCalc_->calc_sh(loc.theta(), loc.phi());
    rscl = _system_->get_aik(I, k) / loc.r();
    chgscl = X[l][h] / loc.r();
    for (n = 0; n < p_; n++)
    {
      for (m = -n; m <= n; m++)
      {
        val = chgscl * _shCalc_->get_result(n, m) + Z(n, m+p_);
        Z.set_val(n, m+p_, val);
      }
      chgscl *= rscl;
    }
  }
  return Z;
}


MyMatrix<cmplx> TMatrix::re_expandX(MyMatrix<cmplx> X,
                                    int I, int k,
                                   int J, int l)
{
  MyMatrix<cmplx> X1, X2, Z;
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE;
  
  X1 = expand_RX(X, I, k, J, l, whichR);
  X2 = expand_SX(X1, I, k, J, l, whichS);
  Z  = expand_RHX(X2, I, k, J, l, whichRH);
  return Z;
}


/*
 re-expand element j of grad(X) with element (i, j) of T. Requires 
 the three components of grad(X)
 */
VecOfMats<cmplx>::type TMatrix::re_expand_gradX(MyMatrix<cmplx> dXdR,
                                                MyMatrix<cmplx> dXdTheta,
                                                MyMatrix<cmplx> dXdPhi,
                                                int I, int k, int J, int l)
{
  MyMatrix<cmplx> x1, x2, z;
  VecOfMats<cmplx>::type Z (3);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE;
  
  // first dA/dR
  x1 = expand_RX(dXdR, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z  = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(0, z);
  
  // dA/dtheta:
  x1 = expand_RX(dXdTheta, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z  = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(1, z);
  
  // dA/dphiL
  x1 = expand_RX(dXdPhi, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z  = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(2, z);
  
  return Z;
}


VecOfMats<cmplx>::type TMatrix::re_expandX_gradT(MyMatrix<cmplx> X,
                                                 int I, int k,
                                                 int J, int l)
{
  MyMatrix<cmplx> x1, x2, z, z1, z2;
  VecOfMats<cmplx>::type Z (3); // output vector
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE;

  
  // first find with respect to dT/dr:
  whichS = DDR;
  x1 = expand_RX(X, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z  = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(0, z);
  
  // dT/dtheta:
  whichR=BASE;
  whichS = BASE;
  whichRH = DDTHETA;
  x1 = expand_RX(X, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z1 = expand_RHX(x2, I, k, J, l, whichRH);
  
  whichRH = BASE;
  whichR = DDTHETA;
  x1 = expand_RX(X, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z2 = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(1, z1 + z2);
  
  // dT/dphi:
  whichR = BASE;
  whichRH = DDPHI;
  x1 = expand_RX(X, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z1 = expand_RHX(x2, I, k, J, l, whichRH);
  
  whichRH = BASE;
  whichR = DDPHI;
  x1 = expand_RX(X, I, k, J, l, whichR);
  x2 = expand_SX(x1, I, k, J, l, whichS);
  z2 = expand_RHX(x2, I, k, J, l, whichRH);
  Z.set_val(2, z1 + z2);
  
  Z = conv_to_cart(Z, I, k, J, l);
  return Z;
}


bool TMatrix::is_Jl_greater(int I, int k, int J, int l)
{
  if (J > I || (I==J && l > k)) return true;
  else return false;
}

// perform first part of T*A and return results
MyMatrix<cmplx> TMatrix::expand_RX(MyMatrix<cmplx> X,
                                   int I, int k, int J, int l,
                                   WhichReEx whichR
                                   )
{

  int n, m, s, map_idx;
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  
  MyMatrix<cmplx> x1(p_, 2*p_ + 1);
  cmplx inter, rval, aval;
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_[map_idx]->isSingular())
      {
        Pt vec = T_[map_idx]->get_TVec();
        if (whichR == BASE)
        {
          aval = X(n, -m+p_); // Not sure if this is right
          if (vec.theta() > M_PI/2.0)
            x1.set_val(n, m+p_, (n%2 == 0 ? aval : -aval));
          else x1.set_val(n, m+p_, X(n, m+p_));  // Not sure if this is right
        } else if (whichR == DDTHETA)
        {
          x1 = expand_dRdtheta_sing(X, I, k, J, l, vec.theta(), false);  // not sure if this is right
          return x1;
        }
        else
        {
          x1 = expand_dRdphi_sing(X, I, k, J, l, vec.theta(), false);  // not sure if this is right
          return x1;
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichR == DDPHI)
            rval = T_[map_idx]->get_dr_dphi_val(n, m, s);
          else if (whichR == DDTHETA)
            rval = T_[map_idx]->get_dr_dtheta_val(n, m, s);
          else
            rval = T_[map_idx]->get_rval(n, m, s);
          
          aval = X(n, s+p_);
          
          inter += rval * aval;
        } // end s
        x1.set_val(n, m+p_, inter);
      }
    } // end m
  } //end n
  return x1;
}

// perform second part of T*A and return results
MyMatrix<cmplx> TMatrix::expand_SX(MyMatrix<cmplx> x1,
                                   int I, int k, int J, int l,
                                   WhichReEx whichS)
{
  
  
  cmplx inter, sval;
  MyMatrix<cmplx> x2(p_, 2*p_ + 1);
  
  int n, m, s, map_idx;
  
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (s = abs(m); s < p_; s++)
      {
        if ( jl_greater )
        {
          if (whichS == DDR) sval = T_[map_idx]->get_dsdr_val(n, s, m);
          else sval = T_[map_idx]->get_sval(n, s, m);
        }
        else
        {
          if (whichS == DDR) sval = T_[map_idx]->get_dsdr_val(s, n, m);
          else sval = T_[map_idx]->get_sval(s, n, m);
        }
        inter += sval * x1(s, m+p_);
      } // end l
      x2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x2;
}
//
// perform third part of T*A and return results
MyMatrix<cmplx> TMatrix::expand_RHX( MyMatrix<cmplx> x2,
                                    int I, int k, int J, int l,
                                    WhichReEx whichRH)
{
  int n, m, s, map_idx;
  cmplx inter, rval;
  MyMatrix<cmplx> z(p_, 2*p_ + 1);
  
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  

  //fill zj:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_[map_idx]->isSingular())
      {
        Pt vec = T_[map_idx]->get_TVec();
        if (whichRH == BASE)
        {
          if (vec.theta() > M_PI/2.0)
            z.set_val(n, m+p_, (n%2 == 0 ? x2(n,-m+p_) : -x2(n,-m+p_)));
          else
            z.set_val(n, m+p_, x2(n,m+p_));
        } else if (whichRH == DDTHETA)
        {
          z = expand_dRdtheta_sing(x2, I, k, J, l, vec.theta(), true);
          return z;
        }
        else
        {
          z = expand_dRdphi_sing(x2, I, k, J, l, vec.theta(), true);
          return z;
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichRH == DDPHI)
            rval = T_[map_idx]->get_dr_dphi_val(n, s, m);
          else if (whichRH == DDTHETA)
            rval = T_[map_idx]->get_dr_dtheta_val(n, s, m);
          else
            rval = T_[map_idx]->get_rval(n, s, m);
          
          inter += conj(rval) * x2(n, s+p_);
        } // end s
        z.set_val(n, m+p_, inter);
      }
    } // end m
  } //end n
  return z;
}

MyMatrix<cmplx> TMatrix::expand_dRdtheta_sing(MyMatrix<cmplx> mat,
                                              int I, int k, int J, int l,
                                              double theta, bool ham)
{
  MyMatrix<cmplx> x(p_, 2*p_ + 1);
  x.set_val( 0, p_, cmplx(0.0, 0.0));
  double rec = (ham ? -1.0 : 1.0);
  int map_idx;
  
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      x.set_val( n, p_, rec*cmplx(2.0*T_[map_idx]->get_prefac_dR_val(n,0,1)*
                                  mat(n, 1+p_).real(), 0.0)); // m = 0
      for (int m = 1; m < n; m++)
      {
        x.set_val( n, m+p_, rec*T_[map_idx]->get_prefac_dR_val(n,m,0)*mat(n,m-1+p_)
                  + rec*T_[map_idx]->get_prefac_dR_val(n,m,1)*mat(n,m+1+p_));
        x.set_val( n,-m+p_, conj( x( n, m+p_)));
      }
      
      x.set_val(n, n+p_,
                rec*T_[map_idx]->get_prefac_dR_val( n, n, 0)*mat(n, n-1+p_));
      x.set_val(n,-n+p_, conj( x( n, n+p_)));
    }
  }
  else
  {
    double s = -1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      x.set_val( n, p_, rec*cmplx(2.0*s*T_[map_idx]->get_prefac_dR_val(n,0,1)*
                                  mat(n,1+p_).real(), 0.0)); // m = 0
      for (int m = 1; m < n; m++)
      {
        x.set_val( n, m+p_, rec*s*
                  (T_[map_idx]->get_prefac_dR_val(n,m,0)*mat(n,-m+1+p_)
                   + T_[map_idx]->get_prefac_dR_val(n,m,1)*mat(n,-m-1+p_)));
        x.set_val( n,-m+p_, conj( x( n, m+p_)));
      }
      
      x.set_val(n, n+p_,rec*s*T_[map_idx]->get_prefac_dR_val(n, n,0)*mat(n,-n+1+p_));
      x.set_val(n,-n+p_, conj( x( n, n+p_)));
    }
  }
  return x;
}

MyMatrix<cmplx> TMatrix::expand_dRdphi_sing(MyMatrix<cmplx> mat,
                                            int I, int k, int J, int l,
                                            double theta, bool ham)
{
  MyMatrix<cmplx> x(p_, 2*p_ + 1);
  x.set_val( 0, p_, cmplx(0.0, 0.0));
  double rec = ((ham && (theta < M_PI/2)) ? -1.0 : 1.0);
  int map_idx;
  
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      x.set_val( n, p_, rec*cmplx(2.0*T_[map_idx]->get_prefac_dR_val(n,0,1)*
                                  mat(n, 1+p_).imag(),0.0));
      for (int m = 1; m < n; m++)
      {
        x.set_val(n, m+p_,
                  rec*(cmplx( 0.0, T_[map_idx]->get_prefac_dR_val(n,m, 0))*mat(n,m-1+p_)
                       - cmplx( 0.0, T_[map_idx]->get_prefac_dR_val(n,m,1))*mat(n,m+1+p_)));
        x.set_val( n, -m+p_, conj(x(n, m+p_)));
      }
      
      x.set_val( n, n+p_,
                rec*cmplx( 0.0, T_[map_idx]->get_prefac_dR_val(n,n,0))*mat(n,n-1+p_));
      x.set_val( n, -n+p_, conj(x( n, n+p_)));
    }
  }
  else
  {
    double s = 1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      x.set_val(n, p_, rec*cmplx(2.0*s*T_[map_idx]->get_prefac_dR_val(n,0,1)*
                                 mat(n,1+p_).imag(),0.0));
      for (int m = 1; m < n; m++)
      {
        x.set_val(n, m+p_,
                  rec*s*(-cmplx(0.0,T_[map_idx]->get_prefac_dR_val(n,m,0))*mat(n,-m+1+p_)
                         + cmplx(0.0, T_[map_idx]->get_prefac_dR_val(n,m,1))*mat(n,-m-1+p_)));
        x.set_val( n, -m+p_, conj(x(n, m+p_)));
      }
      
      x.set_val( n,  n+p_, rec*cmplx(0.0,
                                     -s*T_[map_idx]->get_prefac_dR_val(n,n,0))*mat(n,-n+1+p_));
      x.set_val( n, -n+p_, conj(x(n, n+p_)));
    }
  }
  return x;
}


VecOfMats<cmplx>::type TMatrix::conv_to_cart(VecOfMats<cmplx>::type dZ,
                                             int I, int k, int J, int l)
{
  int m, n, map_idx;
  double the, r, phi;
  VecOfMats<cmplx>::type Zcart (3); vector<double> con1(3), con2(3), con3(3);
  
  bool jl_greater = is_Jl_greater(I, k, J, l);
  if (jl_greater) map_idx = idxMap_[{I, k, J, l}];
  else            map_idx = idxMap_[{J, l, I, k}];
  
  Pt v = T_[map_idx]->get_TVec();
  the = v.theta(); r = v.r(); phi = v.phi();
  
  if (T_[map_idx]->isSingular())
  {
    double cost = (the < M_PI/2) ? 1.0 : -1.0;
    con1 = {  0.0, cost/r,   0.0};
    con2 = {  0.0,    0.0, 1.0/r};
    con3 = { cost,    0.0,   0.0};
  }else
  {
    con1 = {sin(the)*cos(phi),  cos(the)*cos(phi)/r, -sin(phi)/sin(the)/r};
    con2 = {sin(the)*sin(phi),  cos(the)*sin(phi)/r,  cos(phi)/sin(the)/r};
    con3 = {         cos(the),          -sin(the)/r,                  0.0};
  }
  
  for (n = 0; n < 3; n++)
    Zcart[n] = MyMatrix<cmplx> (p_, 2*p_ + 1);
  
  for ( n = 0; n < p_; n++ )
    for ( m = -n; m <= n; m++ )
    {
      (&Zcart[0])->set_val( n, m+p_, con1[0]*dZ[0](n, m+p_)
                           + con1[1]*dZ[1](n, m+p_)+con1[2]*dZ[2](n, m+p_));
      (&Zcart[1])->set_val( n, m+p_, con2[0]*dZ[0](n, m+p_)
                           + con2[1]*dZ[1](n, m+p_)+con2[2]*dZ[2](n, m+p_));
      (&Zcart[2])->set_val( n, m+p_, con3[0]*dZ[0](n, m+p_)
                           + con3[1]*dZ[1](n, m+p_)+con3[2]*dZ[2](n, m+p_));
    }
  return Zcart;
}
//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"

ASolver::ASolver(shared_ptr<BesselCalc> _bcalc,
                  shared_ptr<SHCalc> _shCalc,
                  shared_ptr<System> _sys,
                  shared_ptr<Constants> _consts,
                  const int p,
                  double polz_cutoff)
:p_(p),
N_(_sys->get_n()),
a_avg_(_sys->get_lambda()),
solvedA_(false),
_besselCalc_(_bcalc),
_shCalc_(_shCalc),
polz_cutoff_(polz_cutoff),
_consts_(_consts)
{
  _gamma_ = make_shared<VecOfVecs<double>::type>(N_, MyVector<double> (p_));
  _delta_ = make_shared<VecOfVecs<double>::type>(N_, MyVector<double> (p_));
  _E_ = make_shared<vector<MyExpansion> >(N_, MyExpansion(p_));
  
  _A_ = make_shared<vector<MyExpansion> >(N_, MyExpansion(p_));
  _prevA_ = make_shared<vector<MyExpansion> >(N_,MyExpansion(p_));
  
  _L_ = make_shared<vector<MyExpansion> >(N_, MyExpansion(p_));
  _gradL_ = make_shared<MyVector<MyGradExpansion> >(N_, MyGradExpansion(p_));
  
  _gradT_A_ = make_shared<MyMatrix<MyGradExpansion > > (N_, N_);
  _gradA_ = make_shared<MyMatrix<MyGradExpansion> > (N_, N_);
  _prevGradA_ = make_shared<MyMatrix<MyGradExpansion > > (N_, N_);
  
  _allSh_ = make_shared<vector<vector<MyMatrix<cmplx> > > > (N_,
                                      vector<MyMatrix<cmplx>> (N_));
  
  reset_all(_sys);
}

// perform many iterations of the solution for A
void ASolver::solve_A(double prec, int MAX_POL_ROUNDS)
{
  double scale_dev = (double)(p_*(p_+1)*0.5);
  double cng = scale_dev;
  int ct = 0;
  
  while((cng/scale_dev) > prec)
  {
    iter();
    cng = calc_change();
    if (ct > MAX_POL_ROUNDS) break;
    ct++;
  }
  solvedA_ = true;
  calc_L();
}

void ASolver::solve_gradA(double prec, int MAX_POL_ROUNDS)
{
  assert(solvedA_); // must solve a before this
  double scale_dev = (double)(p_*(p_+1)*0.5*3.0);
  double cng;
  int j, ct;
  
  pre_compute_gradT_A();
  for ( j = 0; j < N_; j++ )
  {
    ct = 0;
    cng = scale_dev;
    while((cng/scale_dev) > prec)
    {
      grad_iter(j);
      cng = calc_grad_change(j);
      
      if (ct > MAX_POL_ROUNDS) break;
      ct++;
    }
  }
  
  calc_gradL();
}


void ASolver::copy_to_prevA()
{
  for (int i=0; i < _A_->size(); i++)
  {
    _prevA_->operator[](i) = _A_->operator[](i);
  }
}

void ASolver::copy_to_prevGradA(int j)
{
  for (int i=0; i < _A_->size(); i++)
  {
      for (int n = 0; n < p_; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          set_prev_dAdr_ni(i, j, n, m, get_dAdr_ni(i, j, n, m));
          set_prev_dAdtheta_ni(i, j, n, m, get_dAdtheta_ni(i, j, n, m));
          set_prev_dAdphi_ni(i, j, n, m, get_dAdphi_ni(i, j, n, m));
        }
      }
  }
}

// one iteration of numerical solution for A (eq 51 in Lotan 2006)
void ASolver::iter()
{
  Pt v;
  int i, j;
  MyExpansion Z, zj(p_), ai;

  copy_to_prevA();
  bool polz(false), interact(false), prev(true);
  for (i = 0; i <  N_; i++)
  {
    interact = false;
    polz = false;
    // relevant re-expansions:
    Z = MyExpansion (p_);
    for (j = 0; j < N_; j++)
    {
      if (i == j) continue;
      
      v = _sys_->get_pbc_dist_vec(i, j);
      if (! _sys_->less_than_cutoff(v) ) continue; // cutoff for interaction

      _sys_->add_J_to_interact_I(i,j);
      interact = true;
      if (v.norm() > polz_cutoff_+_sys_->get_ai(i)+_sys_->get_ai(j)) continue;

      re_expandA(i, j, zj, prev);
      Z += zj;
      _sys_->add_J_to_pol_I(i,j);
      polz = true;
    }
    
    if (polz)
    {
      ai = Z * _delta_->operator[](i);
      ai += _E_->operator[](i);
 
      _A_->operator[](i) = ai * _gamma_->operator[](i);
    } else if (interact)
    {
      _A_->operator[](i) = _E_->operator[](i) * _gamma_->operator[](i);
    }
  }
}

// one iteration of numerical solution for grad(A) WRT j (eq 53 in Lotan 2006)
// Solving for grad_j(A^(i)) by iterating through T^(i,k)
// relevant re-expansions (g prefix means gradient):
void ASolver::grad_iter(int j)
{
  Pt v;
  int i, k;
  MyGradExpansion add(p_);
  
  bool prev(true), polz(false), interact(false); //want to re-expand previous
  copy_to_prevGradA(j);
  for (i = 0; i < N_; i++) // MoleculeAM of interest
  {
    interact = false;
    polz = false;
    MyGradExpansion aij(p_);
    aij = get_gradT_Aij(j, i);
    
    for (k = 0; k < N_; k++) // other MoleculeAMs
    {
      if (k == i) continue;
      v = _sys_->get_pbc_dist_vec(i, k);
      if (! _sys_->less_than_cutoff(v) ) continue;
      interact = true;
      if (v.norm() > polz_cutoff_+_sys_->get_ai(i)+_sys_->get_ai(j)) continue;
      re_expand_gradA(i, k, j, add, prev); // T^(i,k) * grad_j A^(k)
      aij += add;
      polz = true;
    }
   
    if (interact) 
    {
      MyVector<double> gamma_delta = (_gamma_->operator[](i).
                                      mult(_delta_->operator[](i)));
      
      aij.set_dim(0, aij.get_dim(0) * gamma_delta);
      aij.set_dim(1, aij.get_dim(1) * gamma_delta);
      aij.set_dim(2, aij.get_dim(2) * gamma_delta);
      _gradA_->set_val(i, j, aij);
    }
  }
}

double ASolver::calc_grad_change(int wrt)
{
  double change = 0.0;
  int i;
  WhichReEx whichA;
  for (i = 0; i < 3; i++)
  {
    if (i == 0) whichA = DDR;
    else if (i == 1) whichA = DDTHETA;
    else whichA = DDPHI;
    change += calc_change(whichA, wrt);
  }
  return change;
}

double ASolver::calc_change(WhichReEx whichA, int wrt)
{
  int i, k, m;
  double change = 0;
  cmplx prev, curr, a; // intermediate values
  for (i = 0; i < N_; i++)
  {
    for(k = 0; k < p_; k++)
    {
      for(m = 0; m <= k; m++)
      {
        // get value from correct verison of A:
        if (whichA==DDR)
        {
          prev = get_prev_dAdr_ni(i, wrt, k, m);
          curr = get_dAdr_ni(i, wrt, k, m);
        }
        else if (whichA==DDTHETA)
        {
          prev = get_prev_dAdtheta_ni(i, wrt, k, m);
          curr = get_dAdtheta_ni(i, wrt, k, m);
        }
        else if (whichA==DDPHI)
        {
          prev = get_prev_dAdphi_ni(i, wrt, k, m);
          curr = get_dAdphi_ni(i, wrt, k, m);
        }
        else
        {
          prev = get_prevA_ni(i, k, m);
          curr = get_A_ni(i, k, m);
        }
        
        if (prev == cmplx(0.0,0.0) && curr == cmplx(0.0,0.0))
          continue;

        if (fabs(prev.real()) < 1e-30 && fabs(prev.imag()) < 1e-30 )
          a = curr;
        else if (fabs(curr.real()) < 1e-30 && fabs(curr.imag()) < 1e-30)
          a = prev;
        else
          a = 0.5*((prev - curr)/(prev + curr));
        
        change += a.real()*a.real() + a.imag()*a.imag();
      }
    }
  }
  return change;
}

// precompute gradT times A(i,j) for all pairs of MoleculeAMs
void ASolver::pre_compute_gradT_A()
{
  // Solving for grad_j(A^(i)) by iterating through T^(i,k)
  int i, j, k, dim;
  Pt vij, vik;
  double sign;
  
  // relevant re-expansions (g prefix means gradient):
  MyGradExpansion gjT_Ai(p_), gTA(p_);
  bool prev = true; //want to re-expand previous
  for (i = 0; i < N_; i++) // MoleculeAM of interest
  {
    for (j = 0; j < N_; j++) // gradient of interest
    {
      vij = _sys_->get_pbc_dist_vec(i, j);
      gjT_Ai = MyGradExpansion(p_);
      
      if (j == i)
      {
        for (k = 0; k < N_; k++)
        {
          if ( k == i ) continue;
          vik = _sys_->get_pbc_dist_vec(i, k);
          if (! _sys_->less_than_cutoff(vik) ) continue;
          
          re_expandA_gradT( i, k, gTA, prev); // grad_i T^(i,k) A^(k)

          sign = ( k < i ) ? -1.0 : 1.0;
          for (dim = 0; dim < 3; dim++)
            gTA.set_dim(dim, gTA.get_dim(dim)*sign);
        
          gjT_Ai += gTA;
        }
      }
      else if ((j > i) && (_sys_->less_than_cutoff(vij) ))
        re_expandA_gradT( j, i, gjT_Ai, prev); // grad_j T^(j,i) A^(i)
      else if (_sys_->less_than_cutoff(vij) )
      {
        re_expandA_gradT( j, i, gjT_Ai, prev); // grad_j T^(j,i) A^(i)
        for (dim = 0; dim < 3; dim++)
          gjT_Ai.set_dim(dim, gjT_Ai.get_dim(dim)*-1.0);
      }
      _gradT_A_->set_val(i, j, gjT_Ai);
    }
  }
}

//re-expand element j of A with element (i, j) of T and return results
void ASolver::re_expandA(int i, int j, MyExpansion& z, bool prev)
{
  MyExpansion x1(p_), x2(p_);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  expand_RX(  i, j, x1, whichR, whichA, prev);
  expand_SX(  i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z, whichRH);
}


/*
 re-expand element j of grad(A) with element (i, j) of T
 will re-expand the gradient with respect to the wrt input
 */
void ASolver::re_expand_gradA(int i, int j, int wrt,
                              MyGradExpansion & Z, bool prev)
{
//  MyMatrix<cmplx> x1, x2, z;
  MyExpansion x1(p_), x2(p_), z(p_);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=DDR;

  // first dA/dR
  expand_RX(  i, j, x1, whichR, whichA, prev, wrt);
  expand_SX(  i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z, whichRH);
  Z.set_dim(0, z);

  // dA/dtheta:
  whichA = DDTHETA;
  expand_RX(  i, j, x1, whichR, whichA, prev, wrt);
  expand_SX(  i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z, whichRH);
  Z.set_dim(1, z);
  
  // dA/dphiL
  whichA = DDPHI;
  expand_RX(  i, j, x1, whichR, whichA, prev, wrt);
  expand_SX(  i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z, whichRH);
  Z.set_dim(2, z);

}

/*
 Re-expand element j of A with element (i, j) of grad(T) and return
 as a 3-element vector containing the results for the re-expansion
 with each element of grad(T)
 */
void ASolver::re_expandA_gradT(int i, int j, MyGradExpansion &Z, bool prev)
{
//  MyMatrix<cmplx> x1, x2, z, z1, z2;
  MyExpansion x1(p_), x2(p_), z(p_), z1(p_), z2(p_);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  int lowI = i; int hiJ = j;
  if ( i > j ) { lowI = j; hiJ  = i; }
  
  // first find with respect to dT/dr:
  whichS = DDR;
  expand_RX( i, j, x1, whichR, whichA, prev);
  expand_SX( i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z, whichRH);
  Z.set_dim(0, z);
  
  // dT/dtheta:
  whichR=BASE;
  whichS = BASE;
  whichRH = DDTHETA;
  expand_RX( i, j, x1, whichR, whichA, prev);
  expand_SX( i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z1, whichRH);

  whichRH = BASE;
  whichR = DDTHETA;
  expand_RX( i, j, x1, whichR, whichA, prev);
  expand_SX( i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z2, whichRH);
  Z.set_dim(1, z1 + z2);
  
  // dT/dphi:
  whichR = BASE;
  whichRH = DDPHI;
  expand_RX( i, j, x1, whichR, whichA, prev);
  expand_SX( i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z1, whichRH);
  
  whichRH = BASE;
  whichR = DDPHI;
  expand_RX( i, j, x1, whichR, whichA, prev);
  expand_SX( i, j, x1, x2, whichS);
  expand_RHX( i, j, x2, z2, whichRH);
  Z.set_dim(2, z1 + z2);

  Z = conv_to_cart(Z, i, j);
}


// perform first part of T*A and return results
void ASolver::expand_RX(int i, int j, MyExpansion& x1, WhichReEx whichR,
                                   WhichReEx whichA, bool prev, int wrt)
{
  int n, m, s, lowI, hiJ;
  cmplx inter, rval, aval;
  
  lowI = i; hiJ = j;
  if ( i > j ) { lowI = j; hiJ  = i; }
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = 0; m <= n; m++)
    {
      inter  = 0;
      if (T_(lowI,hiJ).isSingular())
      {
        Pt vec = T_(lowI, hiJ).get_TVec();
        if (whichR == BASE)
        {
          aval = which_aval(whichA, prev, j, n, -m, wrt);
          if (vec.theta() > M_PI/2.0)
            x1.set_val_cmplx(n, m, (n%2 == 0 ? aval : -aval));
          else x1.set_val_cmplx(n, m, which_aval(whichA, prev, j, n, m, wrt));
        } else if (whichR == DDTHETA)
        {
          expand_dRdtheta_sing(i, j, vec.theta(), x1, false);
        }
        else
        {
          expand_dRdphi_sing(i, j, vec.theta(), x1, false);
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichR == DDPHI)
            rval = T_(lowI, hiJ).get_dr_dphi_val(n, m, s);
          else if (whichR == DDTHETA)
            rval = T_(lowI, hiJ).get_dr_dtheta_val(n, m, s);
          else
            rval = T_(lowI, hiJ).get_rval(n, m, s);
          
          aval = which_aval(whichA, prev, j, n, s, wrt);
          inter += rval * aval;
        } // end s
        x1.set_val_cmplx(n, m, inter);
      }
    } // end m
  } //end n
}

// perform second part of T*A and return results
void ASolver::expand_SX(int i, int j, MyExpansion x1, MyExpansion& x2,
                        WhichReEx whichS)
{
  int n, m, l, lowI, hiJ;
  cmplx inter, sval;
  
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = 0; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        if ( i < j )
        {
          if (whichS == DDR) sval = T_(lowI, hiJ).get_dsdr_val(n, l, m);
          else sval = T_(lowI, hiJ).get_sval(n, l, m);
        }
        else
        {
          if (whichS == DDR) sval = T_(lowI, hiJ).get_dsdr_val(l, n, m);
          else sval = T_(lowI, hiJ).get_sval(l, n, m);
        }
        inter += sval * x1.get_cmplx(l, m);
      } // end l
      x2.set_val_cmplx(n, m, inter);
    } // end m
  } //end n
}


// perform third part of T*A and return results
void ASolver::expand_RHX(int i, int j, MyExpansion x2, MyExpansion &z,
                         WhichReEx whichRH)
{
  int n, m, s, lowI, hiJ;
  cmplx inter, rval;
  double mult;
  
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  
  //fill zj:
  for (n = 0; n < p_; n++)
  {
    for (m = 0; m <= n; m++)
    {
      inter  = 0;
      if (T_(lowI,hiJ).isSingular())
      {
        Pt vec = T_(lowI, hiJ).get_TVec();
        if (whichRH == BASE)
        {
          if (vec.theta() > M_PI/2.0)
            z.set_val_cmplx(n, m, (n%2 == 0 ? x2.get_cmplx(n,-m) :
                                   -x2.get_cmplx(n,-m)));
          else
            z.set_val_cmplx(n, m, x2.get_cmplx(n,m));
        } else if (whichRH == DDTHETA)
        {
          expand_dRdtheta_sing(lowI, hiJ, vec.theta(), x2, z, true);
        }
        else
        {
          expand_dRdphi_sing(lowI, hiJ, vec.theta(), x2, z, true);
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichRH == DDPHI)
            rval = T_(lowI, hiJ).get_dr_dphi_val(n, s, m);
          else if (whichRH == DDTHETA)
            rval = T_(lowI, hiJ).get_dr_dtheta_val(n, s, m);
          else
            rval = T_(lowI, hiJ).get_rval(n, s, m);
          
          mult = (s==0) ? 1.0 : 2.0;
          inter += conj(rval) * x2.get_cmplx(n, s);
        } // end s
        z.set_val_cmplx(n, m, inter);
      }
    } // end m
  } //end n
}

void ASolver::expand_dRdtheta_sing(int i, int j, double theta,
                                   MyExpansion& mat, bool ham)
{
  int lowI, hiJ;
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  expand_dRdtheta_sing(lowI, hiJ, theta, _prevA_->operator[](j), mat, ham);
}

void ASolver::expand_dRdtheta_sing(int i, int j, double theta,
                                   MyExpansion mat_in,
                                   MyExpansion& mat_out, bool ham)
{
  mat_out.set_val_cmplx( 0, 0, cmplx(0.0, 0.0));
  double rec = (ham ? -1.0 : 1.0);
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      mat_out.set_val_cmplx(n, 0, rec*
                            cmplx(2.0*T_(i,j).get_prefac_dR_val(n,0,1)*
                            mat_in(n, 1, REAL), 0.0)); // m = 0
      
      for (int m = 1; m < n; m++)
      {
        mat_out.set_val_cmplx( n, m, rec*T_(i,j).get_prefac_dR_val(n,m,0)*
                              mat_in.get_cmplx(n,m-1) + rec*T_(i,j).get_prefac_dR_val(n,m,1)*mat_in.get_cmplx(n,m+1));
      }

      mat_out.set_val_cmplx(n, n, rec*T_(i,j).get_prefac_dR_val( n, n, 0)
                            *mat_in.get_cmplx(n, n-1));
    }
  }
  else
  {
    double s = -1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      mat_out.set_val_cmplx( n, 0, rec*cmplx(2.0*s*T_(i,j).get_prefac_dR_val(n,0,1)*
                                  mat_in.get_cmplx(n,1).real(), 0.0)); // m = 0
      for (int m = 1; m < n; m++)
      {
        mat_out.set_val_cmplx( n, m, rec*s*
                  (T_(i,j).get_prefac_dR_val(n,m,0)*mat_in.get_cmplx(n,-m+1)
                   + T_(i,j).get_prefac_dR_val(n,m,1)*mat_in.get_cmplx(n,-m-1)));
      }
      
      mat_out.set_val_cmplx(n, n,rec*s*T_(i,j).get_prefac_dR_val(n, n,0)
                            *mat_in.get_cmplx(n,-n+1));
    }
  }
}

void ASolver::expand_dRdphi_sing(int i, int j, double theta,
                                 MyExpansion& mout, bool ham)
{
  int lowI, hiJ;
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  expand_dRdphi_sing(lowI, hiJ, theta, _prevA_->operator[](j), mout, ham);
}

void ASolver::expand_dRdphi_sing(int i, int j, double theta, MyExpansion mat_in,
                                 MyExpansion&mat_out, bool ham)
{
  mat_out.set_val_cmplx( 0, 0, cmplx(0.0, 0.0));
  double rec = ((ham && (theta < M_PI/2)) ? -1.0 : 1.0);
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      mat_out.set_val_cmplx( n, 0, rec*cmplx(2.0*T_(i,j).get_prefac_dR_val(n,0,1)*
                                  mat_in(n, 1, IMAG),0.0));
      for (int m = 1; m < n; m++)
      {
        mat_out.set_val_cmplx(n, m, rec*(cmplx(0.0, T_(i,j).get_prefac_dR_val(n,m, 0))*mat_in.get_cmplx(n,m-1) - cmplx( 0.0, T_(i,j).get_prefac_dR_val(n,m,1))*mat_in.get_cmplx(n,m+1)));
      }
      
      mat_out.set_val_cmplx( n, n, rec*cmplx( 0.0, T_(i,j).
                                             get_prefac_dR_val(n,n,0))
                            *mat_in.get_cmplx(n,n-1));
    }
  }
  else
  {
    double s = 1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      mat_out.set_val_cmplx(n, 0, rec*cmplx(2.0*s*T_(i,j).get_prefac_dR_val(n,0,1)*
                                 mat_in.get_cmplx(n,1).imag(),0.0));
      for (int m = 1; m < n; m++)
      {
        mat_out.set_val_cmplx(n, m,rec*s*(-cmplx(0.0,T_(i,j).
                                                 get_prefac_dR_val(n,m,0))
                                          *mat_in.get_cmplx(n,-m+1)
                                          + cmplx(0.0, T_(i,j).
                                                  get_prefac_dR_val(n,m,1))
                                          *mat_in.get_cmplx(n,-m-1)));
      }
      
      mat_out.set_val_cmplx( n,  n, rec*cmplx(0.0, -s*T_(i,j).
                                              get_prefac_dR_val(n,n,0))
                            *mat_in.get_cmplx(n,-n+1));
    }
  }
}

/*
 Return correct aval, given A enum and prev flag, mol index, and nm index
 */
cmplx ASolver::which_aval(WhichReEx whichA, bool prev, int i, int n,
                          int m, int wrt)
{
  cmplx aval;
  if (prev)
  {
    if (whichA == DDR)
      aval = get_prev_dAdr_ni(i, wrt, n, m);
    else if (whichA == DDTHETA)
      aval = get_prev_dAdtheta_ni(i, wrt, n, m);
    else if (whichA == DDPHI)
      aval = get_prev_dAdphi_ni(i, wrt, n, m);
    else
      aval = get_prevA_ni(i, n, m);
  }
  else
  {
    if (whichA == DDR)
      aval = get_dAdr_ni(i, wrt, n, m);
    else if (whichA == DDTHETA)
      aval = get_dAdtheta_ni(i, wrt, n, m);
    else if (whichA == DDPHI)
      aval = get_dAdphi_ni(i, wrt, n, m);
    else
      aval = get_A_ni(i, n, m);
  }
  return aval;
}

/*
 Precompute spherical harmonics for E calculation
 */
void ASolver::pre_compute_all_sh()
{
  int i;
  for (i = 0; i < N_; i++)
    (*_allSh_)[i] = calc_mol_sh(_sys_->get_moli(i));
}

/*
 Calculate the SH matrix for every charge in a MoleculeAM
 */
vector<MyMatrix<cmplx> > ASolver::calc_mol_sh(MoleculeAM mol)
{
  vector<MyMatrix<cmplx> > vout;
  vout.reserve(mol.get_m());
  int j;
  double theta, phi;
  Pt pt;
  for (j = 0; j < mol.get_m(); j++)
  {
    pt = mol.get_posj(j);
    theta = pt.theta();
    phi = pt.phi();
    _shCalc_->calc_sh(theta, phi);
    vout.push_back(_shCalc_->get_full_result());
  }
  return vout;
}

/*
 Equation 19--Lotan 2006 page 544
 */
double ASolver::calc_indi_gamma(int i, int n)
{
  double g;  // output
  double kap = _consts_->get_kappa();
  double ai = _sys_->get_ai(i);
  double eps_p = _consts_->get_dielectric_prot();
  double eps_s = _consts_->get_dielectric_water();
  // all bessel function k
  vector<double> bk_all = _besselCalc_->calc_mbfK(n+2, kap*ai);
  double bk1 = bk_all[n];  // bessel k at n
  double bk2 = bk_all[n+1];   // bessel k at n+1
  
  double n_dub = (double) n;
  g = (2.0*n_dub + 1.0) * exp(kap*ai);
  g = g / (((2.0*n_dub + 1.0)*bk2) + (n_dub*bk1*((eps_p / eps_s) - 1.0)));
  return g;
}

/*
 Equation 20--Lotan 2006 page 544
 */
double ASolver::calc_indi_delta(int i, int n)
{
  double d;  // output
  double kap = _consts_->get_kappa();
  double ai = _sys_->get_ai(i);  // radius
  double eps_p = _consts_->get_dielectric_prot();
  double eps_s = _consts_->get_dielectric_water();
  // all bessel function I:
  vector<double> bi_all = _besselCalc_->calc_mbfI(n+2, kap*ai);
  double bi1 = bi_all[n];  // bessel i at n
  double bi2 = bi_all[n+1];   // bessel i at n+1
  
  double n_dub = (double) n;
  d = (kap*kap*ai*ai) * (bi2 / (2.0*n_dub + 3.0));
  d += (n_dub * bi1 * (1.0 - (eps_p / eps_s)));
  d *= (ai * pow(ai/a_avg_, 2.0*n_dub)) / (2.0*n_dub+1.0);
  return d;
}

/*
 Calculates an E^(i)_(n,m) value
 Equation 13--Lotan 2006 page 543
 */
cmplx ASolver::calc_indi_e(int i, int n, int m)
{
  cmplx e = 0.0;
  int j;
  double q, rho, lambda;
  for (j = 0; j < _sys_->get_Mi(i); j++)
  {
    q = _sys_->get_qij(i, j);
    rho = _sys_->get_posij(i, j).r();
    lambda = pow(_sys_->get_lambda(), n);
    // q_ij * (rho_ij)^n * Y_(n,m)(theta_ij, phi_ij):
    cmplx all_sh_acc = (*_allSh_)[i][j](n, abs(m));
    if ( m < 0 )
      all_sh_acc = conj( all_sh_acc );

    e += q * ( pow( rho, n ) / lambda ) * all_sh_acc;
  }
  return e;
}


/*
 Constructs the gamma matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_gamma()
 */
void ASolver::compute_gamma()
{
  int i, j;
  for (i = 0; i < N_; i++)
  {
    for(j = 0; j < p_; j++)
      _gamma_->operator[](i).set_val(j, calc_indi_gamma(i, j));
  }
}


/*
 Constructs the delta matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_delta()
 */
void ASolver::compute_delta()
{
  int i, j;
  
  for (i = 0; i < N_; i++)
  {
    for(j = 0; j < p_; j++)
      _delta_->operator[](i).set_val(j, calc_indi_delta(i, j));
  }
}

/*
 Constructs the E vector, which contains a matrix for each MoleculeAM
 that defines the multipole expansion of that MoleculeAM. The values
 of the inner matrices are calculated in calc_indi_e()
 */
void ASolver::compute_E()
{
  int i, n, m;
  for (i = 0; i < N_; i++)
  {
    // m goes from -n to n so you need 2*p columns:
    for (n = 0; n < p_; n++)
    {
      for (m = 0; m <= n; m++)
        _E_->operator[](i).set_val_cmplx(n, m, calc_indi_e(i, n, m));
    }
  }
}

/*
  Calculate the T matrix (i.e the re expansion coefficients along
  every inter-molecular vector (distance between every molecular center
 */
void ASolver::compute_T()
{
  int i, j;
  Pt v, ci, cj;  // inter molecular vector
  
  for (i = 0; i < N_; i++)
  {
    for (j = i+1; j < N_; j++)
    {
      if (i == j) continue; //PBC
      v = _sys_-> get_pbc_dist_vec(i, j);
      if (! _sys_->less_than_cutoff(v)) continue;
       
      // calculate spherical harmonics for inter molecular vector:
      double kappa = _consts_->get_kappa();
      _shCalc_->calc_sh(v.theta(), v.phi());
      vector<double> besselK = _besselCalc_->calc_mbfK(2*p_, kappa * v.r());
      T_.set_val(i, j, ReExpCoeffs(p_, v, _shCalc_->get_full_result(),
                                   besselK, _reExpConsts_,
                                   {kappa,kappa}, {_sys_->get_lambda()},true));
    }
  }
}

// Initialize A to Gamma * E
void ASolver::init_A()
{
  int i;
  MyMatrix<cmplx> ai;
  for (i = 0; i < N_; i++)
  {
    _A_->operator[](i) = _E_->operator[](i) * _gamma_->operator[](i);
  }
}

// Initialize grad(A) matrix to the zero matrix
void ASolver::init_gradA()
{
  int i, j; //, k;
//  VecOfMats<cmplx>::type ga_ij;
//  MyMatrix<cmplx> wrt;
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      _gradA_->set_val(i, j, MyGradExpansion(p_));
    }
  }
}


/*
 Calculate L, equation (16) of Lotan 2006
 */
void ASolver::calc_L()
{
  int i, j;
  Pt v;
  MyExpansion expand(p_);
  for (i = 0; i < N_; i++)
  {
    _L_->operator[](i) = MyExpansion (p_);
    for (j = 0; j < N_; j++)
    {
      if (j == i) continue;
      v = _sys_->get_pbc_dist_vec(i, j);
      if (! _sys_->less_than_cutoff(v) ) continue;
      
      re_expandA(i, j, expand);
      _L_->operator[](i) += expand;
    }
  }
}

void ASolver::calc_gradL()
{
  int i, k;
  Pt v;
  MyGradExpansion inner1(p_), inner2(p_);
  
  for (i = 0; i < N_; i++) // MoleculeAM of interest
  {
    inner1 = get_gradT_Aij( i, i);
    for (k = 0; k < N_; k++) // other MoleculeAMs
    {
      if (k == i) continue;
      v = _sys_->get_pbc_dist_vec(i, k);
      if (! _sys_->less_than_cutoff(v) ) continue;
      
      re_expand_gradA(i, k, i, inner2, false); // T^(i,k) * grad_j A^(k)
      inner1 += inner2;
    }
    _gradL_->operator[](i) = inner1;
  }
}

/*
 Convert derivatives from spherical to cartesian coords
 */
MyGradExpansion ASolver::conv_to_cart(MyGradExpansion& dZ, int i, int j)
{
  int lowI(i), hiJ(j);
  double the, r, phi;
  MyGradExpansion Zcart(p_);
  vector<double> con1(3), con2(3), con3(3);

  if ( i > j ) { lowI = j; hiJ  = i; }
  
  Pt v = T_(lowI, hiJ).get_TVec();
  the = v.theta(); r = v.r(); phi = v.phi();

  if (T_(lowI, hiJ).isSingular())
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
  
  Zcart.set_dim(0, dZ.get_dim(0)*con1[0] + dZ.get_dim(1)*con1[1] +
                dZ.get_dim(2)*con1[2]);
  Zcart.set_dim(1, dZ.get_dim(0)*con2[0] + dZ.get_dim(1)*con2[1] +
                dZ.get_dim(2)*con2[2]);
  Zcart.set_dim(2, dZ.get_dim(0)*con3[0] + dZ.get_dim(1)*con3[1] +
                dZ.get_dim(2)*con3[2]);
  return Zcart;
}

/*
 Print function for E
 */
void ASolver::print_Ei( int i, int p)
{
  cout << "This is my E for MoleculeAM " << i << " poles " << p <<  endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_E_ni( i, n, m).real();
      double im = get_E_ni( i, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9)  << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for A
 */
void ASolver::print_Ai( int i, int p)
{
  cout << "This is my A for MoleculeAM " << i << " poles " << p <<  endl;
  for (int n = 0; n < 5; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_A_ni( i, n, m).real();
      double im = get_A_ni( i, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << " (" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dx of MoleculeAM i wrt j
 */
void ASolver::print_dAidx( int i, int j, int p)
{
  cout << "dAdx for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdx_ni(i, j, n, m).real();
      double im = get_dAdx_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dy of MoleculeAM i wrt j
 */
void ASolver::print_dAidy( int i, int j, int p)
{
  cout << "dAdy for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdy_ni(i, j, n, m).real();
      double im = get_dAdy_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dz of MoleculeAM i wrt j
 */
void ASolver::print_dAidz( int i, int j, int p)
{
  cout << "dAdz for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdz_ni(i, j, n, m).real();
      double im = get_dAdz_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA of MoleculeAM i wrt j
 */
void ASolver::print_dAi( int i, int j, int p)
{
  print_dAidx( i, j, p);
  print_dAidy( i, j, p);
  print_dAidz( i, j, p);
}


void ASolver::reset_all(shared_ptr<System> _sys)
{
  _sys_ = _sys;
  T_  = MyMatrix<ReExpCoeffs>(_sys->get_n(), _sys->get_n());
  _reExpConsts_ = make_shared<ReExpCoeffsConstants>(_consts_->get_kappa(),
                                                    _sys_->get_lambda(), p_);
  int n, n1;
  for ( n = 0; n < N_; n++ )
  {
    _gamma_->operator[](n) = MyVector<double> (p_);
    _delta_->operator[](n) = MyVector<double> (p_);
    _E_->operator[](n) = MyExpansion (p_);
    
    _A_->operator[](n) = MyExpansion (p_);
    _prevA_->operator[](n) = MyExpansion (p_);
    
    _allSh_->operator[](n) = vector<MyMatrix<cmplx>>
                              (N_, MyMatrix<cmplx> (2*p_, 2*p_));
    
    for ( n1 = 0; n1 < N_; n1++ )
    {
      _gradT_A_->operator()(n, n1) = MyGradExpansion(p_);
      _gradA_->operator()(n, n1) = MyGradExpansion(p_);
      _prevGradA_->operator()(n, n1) = MyGradExpansion(p_);
    }
  }
  
  solvedA_ = false;
  
  // precompute all SH:
  pre_compute_all_sh();
  
  //precompute gamma, delta, E and T:
  compute_T();
  compute_gamma();
  compute_delta();
  compute_E();
  
  init_A();
  init_gradA();
}

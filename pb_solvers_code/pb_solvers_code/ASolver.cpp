//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"
#include <iostream>

ASolver::ASolver(const int N, const int p, BesselCalc bcalc,
                 SHCalc shCalc, System sys)
:p_(p), a_avg_(sys.get_lambda()), gamma_(N),delta_(N), E_(N), N_(sys.get_n()),
T_ (N_, N_), A_(N), gradA_(N, N), prevA_(N), prevGradA_(N, N),
reExpConsts_(sys.get_consts().get_kappa(), sys.get_lambda(), p)
{
  _sys_ = make_shared<System> (sys);
  _besselCalc_ = make_shared<BesselCalc> (bcalc);
  _shCalc_ = make_shared<SHCalc> (shCalc);
  
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

// perform many iterations of the solution for A
void ASolver::solve_A(double prec)
{
  int MAX_POL_ROUNDS = 500;
  prevA_ = A_;
  iter();
  
  double scale_dev = (double)(p_*(p_+1)*0.5);
  double cng = scale_dev;
  int ct = 0;
  
  while((cng/scale_dev) > prec)
  {
    iter();
    cng = calc_change();
    if (ct > MAX_POL_ROUNDS*N_)
    {
      cout << "Polz doesn't converge! dev="<< cng << " " << ct << endl;
      exit(0);
    }
    ct++;
  }
}

// one iteration of numerical solution for A (eq 51 in Lotan 2006)
void ASolver::iter()
{
  
  int i, j;
  MyMatrix<cmplx> Z, zj, ai;
  prevA_ = A_;
  
  for (i = 0; i <  N_; i++)
  {
    // relevant re-expansions:
    Z = MyMatrix<cmplx> (p_, 2*p_ + 1);
    for (j = 0; j < N_; j++)
    {
      if (i == j) continue;
      zj = re_expandA(i, j);
      Z += zj;
    }

    ai = delta_[i] * Z;
    ai += E_[i];
    ai = gamma_[i] * ai;
    A_.set_val(i, ai);
  }
}


// one iteration of numerical solution for grad(A) (eq 53 in Lotan 2006)
void ASolver::grad_iter()
{
  // Solving for grad_j(A^(i)) by iterating through T^(i,k)
  int i, j, k;
  prevGradA_ = gradA_;
  
  // relevant re-expansions (g prefix means gradient):
  VecOfMats<cmplx>::type gTA, TgA;
  
  VecOfMats<cmplx>::type aij;
  VecOfMats<cmplx>::type add;
  MyMatrix<cmplx> gamma_delta;
  
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      aij = VecOfMats<cmplx>::type (3);
      for (k = 0; k < N_; k++)
      {
        if (k == i) continue;
        gTA = re_expandA_gradT(i, k);
        TgA = re_expand_gradA(i, k, j);
        add = gTA + TgA;
        aij += add;
      }
      gamma_delta = gamma_[i] * delta_[i];
      aij.set_val(0, aij[0] * gamma_delta);
      aij.set_val(1, aij[1] * gamma_delta);
      aij.set_val(2, aij[2] * gamma_delta);
      gradA_.set_val(i, j, aij);
    }
  }
  
}


double ASolver::calc_change()
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
        prev = get_prevA_ni(i, k, m);
        curr = get_A_ni(i, k, m);
        
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

void ASolver::pre_compute_all_sh()
{
  all_sh.reserve(N_);
  int i;
  
  for (i = 0; i < N_; i++)
  {
    all_sh.push_back(calc_mol_sh(_sys_->get_molecule(i)));
  }
  
}

//re-expand element j of A with element (i, j) of T and return results
MyMatrix<cmplx> ASolver::re_expandA(int i, int j)
{
  MyMatrix<cmplx> x1, x2, z;
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z  = expand_RHX( i, j, x2, whichRH);
  return z;
}

// re-expand element j of grad(A) with element (i, j) of T
// will re-expand the gradient with respect to the wrt input
VecOfMats<cmplx>::type ASolver::re_expand_gradA(int i, int j, int wrt)
{
  MyMatrix<cmplx> x1, x2, z;
  VecOfMats<cmplx>::type Z (3);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=DDR;
  
  // first dA/dR
  x1 = expand_RX(  i, j, whichR, whichA, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z = expand_RHX( i, j, x2, whichRH);
  Z.set_val(0, z);
  
  // dA/dtheta:
  whichA = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z = expand_RHX( i, j, x2, whichRH);
  Z.set_val(1, z);
  
  // dA/dphiL
  whichA = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z = expand_RHX( i, j, x2, whichRH);
  Z.set_val(2, z);
  
  return Z;
}

/*
 Re-expand element j of A with element (i, j) of grad(T) and return
 as a 3-element vector containing the results for the re-expansion
 with each element of grad(T)
 */
VecOfMats<cmplx>::type ASolver::re_expandA_gradT(int i, int j)
{
  MyMatrix<cmplx> x1, x2, z, z1, z2;
  VecOfMats<cmplx>::type Z (3); // output vector
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  // first find with respect to dT/dr:
  whichS = DDR;
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z = expand_RHX( i, j, x2, whichRH);
  Z.set_val(0, z);
  
  // dT/dtheta:
  whichS = BASE;
  whichRH = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z1  = expand_RHX( i, j, x2, whichRH);
  whichRH = BASE;
  whichR = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z2  = expand_RHX( i, j, x2, whichRH);
  Z.set_val(1, z1 + z2);
  
  // dT/dphi:
  whichR = BASE;
  whichRH = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z1  = expand_RHX( i, j, x2, whichRH);
  whichRH = BASE;
  whichR = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA);
  x2 = expand_SX(  i, j, x1, whichS);
  z2  = expand_RHX( i, j, x2, whichRH);
  Z.set_val(2, z1 + z2);
  
  return Z;
}


// perform first part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RX(int i, int j, WhichReEx whichR,
                                   WhichReEx whichA, int wrt)
{
  int n, m, s;
  cmplx inter;
  MyMatrix<cmplx> x1(p_, 2*p_ + 1);
  cmplx rval;
  cmplx aval;
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_(i,j).isSingular())
      {
        Pt vec = T_(i,j).get_TVec();
        if (whichA == DDR)
          aval = get_dAdr_ni(j, wrt, n, -m);
        else if (whichA == DDTHETA)
          aval = get_dAdtheta_ni(j, wrt, n, -m);
        else if (whichA == DDPHI)
          aval = get_dAdphi_ni(j, wrt, n, -m);
        else
          aval = get_A_ni(j, n, -m);
        
        
        if (vec.theta() > M_PI/2.0)
          inter = (n % 2 == 0 ? aval : -aval);
        else
          inter = aval;
        
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichR == DDPHI)
            rval = T_(i, j).get_dr_dphi_val(n, m, s);
          else if (whichR == DDTHETA)
            rval = T_(i, j).get_dr_dtheta_val(n, m, s);
          else
            rval = T_(i, j).get_rval(n, m, s);
          
          if (whichA == DDR)
            aval = get_dAdr_ni(j, wrt, n, s);
          else if (whichA == DDTHETA)
            aval = get_dAdtheta_ni(j, wrt, n, s);
          else if (whichA == DDPHI)
            aval = get_dAdphi_ni(j, wrt, n, s);
          else
            aval = get_A_ni(j, n, s);
          
          inter += rval * aval;
        } // end s
      }
      x1.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x1;
}

// perform second part of T*A and return results
MyMatrix<cmplx> ASolver::expand_SX(int i, int j, MyMatrix<cmplx> x1,
                                   WhichReEx whichS)
{
  int n, m, l;
  cmplx inter;
  MyMatrix<cmplx> x2(p_, 2*p_ + 1);
  cmplx sval;
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = 0; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        if (whichS == DDR) sval = T_(i, j).get_dsdr_val(n, l, m);
        else sval = T_(i, j).get_sval(n, l, m);
        inter += sval * x1(l, m+p_);
      } // end l
      x2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x2;
}

// perform third part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RHX(int i, int j, MyMatrix<cmplx> x2,
                                    WhichReEx whichRH)
{
  int n, m, s;
  cmplx inter;
  MyMatrix<cmplx> z(p_, 2*p_ + 1);
  cmplx rval;
  
  //fill zj:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (s = -n; s <= n; s++)
      {
        if (whichRH == DDPHI)
          rval = T_(i, j).get_dr_dphi_val(n, m, s);
        else if (whichRH == DDTHETA)
          rval = T_(i, j).get_dr_dtheta_val(n, m, s);
        else
          rval = T_(i, j).get_rval(n, m, s);
        inter += conj(rval) * x2(n, s+p_);
      } // end s
      z.set_val(n, m+p_, inter);
    } // end m
  } //end n
  
  return z;
}

/*
 Calculate the SH matrix for every charge in a molecule
 */
vector<MyMatrix<cmplx> > ASolver::calc_mol_sh(Molecule mol)
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
cmplx ASolver::calc_indi_gamma(int i, int n)
{
  double g;  // output
  double kap = _sys_->get_consts().get_kappa();
  double ai = _sys_->get_ai(i);
  double eps_p = _sys_->get_consts().get_dielectric_prot();
  double eps_s = _sys_->get_consts().get_dielectric_water();
  // all bessel function k
  vector<double> bk_all = _besselCalc_->calc_mbfK(n+2, kap*ai);
  double bk1 = bk_all[n];  // bessel k at n
  double bk2 = bk_all[n+1];   // bessel k at n+1
  
  double n_dub = (double) n;
  g = (2.0*n_dub + 1.0) * exp(kap*ai);
  g = g / (((2.0*n_dub + 1.0)*bk2) + (n_dub*bk1*((eps_p / eps_s) - 1.0)));
  return cmplx(g, 0);
}

/*
 Equation 20--Lotan 2006 page 544
 */
cmplx ASolver::calc_indi_delta(int i, int n)
{
  double d;  // output
  double kap = _sys_->get_consts().get_kappa();
  double ai = _sys_->get_ai(i);  // radius
  double eps_p = _sys_->get_consts().get_dielectric_prot();
  double eps_s = _sys_->get_consts().get_dielectric_water();
  // all bessel function I:
  vector<double> bi_all = _besselCalc_->calc_mbfI(n+2, kap*ai);
  double bi1 = bi_all[n];  // bessel i at n
  double bi2 = bi_all[n+1];   // bessel i at n+1
  
  double n_dub = (double) n;
  d = (kap*kap*ai*ai) * (bi2 / (2.0*n_dub + 3.0));
  d += (n_dub * bi1 * (1.0 - (eps_p / eps_s)));
  d *= (ai * pow(ai/a_avg_, 2.0*n_dub)) / (2.0*n_dub+1.0);
  return cmplx(d, 0);
}

/*
 Calculates an E^(i)_(n,m) value
 Equation 13--Lotan 2006 page 543
 */
cmplx ASolver::calc_indi_e(int i, int n, int m)
{
  cmplx e = 0.0;
  int j;
  double q, rho;
  for (j = 0; j < _sys_->get_Mi(i); j++)
  {
    q = _sys_->get_qij(i, j);
    rho = _sys_->get_posij(i, j).r();
    // q_ij * (rho_ij)^n * Y_(n,m)(theta_ij, phi_ij):
    cmplx all_sh_acc = all_sh[i][j](n, abs(m));
    if ( m < 0 )
      all_sh_acc = conj( all_sh_acc );

    e += q * pow( rho, n ) * all_sh_acc;
    
  }
  return e;
}


/*
 Constructs the gamma matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_gamma()
 */
void ASolver::compute_gamma()
{
  MyMatrix<cmplx> gi;
  int i, j;
  for (i = 0; i < N_; i++)
  {
    gi = MyMatrix<cmplx> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      gi.set_val(j, j, calc_indi_gamma(i, j));
    }
    gamma_.set_val(i, gi);
  }
}


/*
 Constructs the delta matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_delta()
 */
void ASolver::compute_delta()
{
  MyMatrix<cmplx> di;
  int i, j;
  
  for (i = 0; i < N_; i++)
  {
    di = MyMatrix<cmplx> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      di.set_val(j, j, calc_indi_delta(i, j));
    }
    delta_.set_val(i, di);
  }
}

/*
 Constructs the E vector, which contains a matrix for each molecule
 that defines the multipole expansion of that molecule. The values
 of the inner matrices are calculated in calc_indi_e()
 */
void ASolver::compute_E()
{
  int i, n, m;
  MyMatrix<cmplx> ei;
  for (i = 0; i < N_; i++)
  {
    // m goes from -n to n so you need 2*p columns:
    ei = MyMatrix<cmplx>(p_, 2*p_ + 1);
    for (n = 0; n < p_; n++)
    {
      for (m = -n; m <= n; m++)
      {
        ei.set_val(n, m + p_, calc_indi_e(i, n, m));
      }
    }
    E_.set_val(i, ei);
  }
}

void ASolver::compute_T()
{
  int i, j;
  // Calculate the T matrix (i.e the re expansion coefficients along
  // every inter-molecular vector (distance between every molecular center
  Pt v, ci, cj;  // inter molecular vector
  
  
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      if (i == j) continue;
      ci = _sys_->get_centeri(i);
      cj = _sys_->get_centeri(j);
      v = ci - cj;
      // calculate spherical harmonics for inter molecular vector:
      double kappa = _sys_->get_consts().get_kappa();
      _shCalc_->calc_sh(v.theta(), v.phi());
      vector<double> besselK = _besselCalc_->calc_mbfK(2*p_, kappa * v.r());
      T_.set_val(i, j, ReExpCoeffs(p_, v, _shCalc_->get_full_result(),
                                   besselK, reExpConsts_,
                                   kappa, _sys_->get_lambda()));
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
    // m goes from -n to n so you need 2*p columns:
    ai = MyMatrix<cmplx>(p_, 2*p_ + 1);
    ai = gamma_[i] * E_[i];
    
    A_.set_val(i, ai);
  }
  
}

//inialize grad(A) matrix to the zero matrix
void ASolver::init_gradA()
{
  int i, j, k;
  VecOfMats<cmplx>::type ga_ij;
  MyMatrix<cmplx> wrt;
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      ga_ij = VecOfMats<cmplx>::type (3);
      for (k = 0; k < 3; k++)
      {
        wrt = MyMatrix<cmplx>(p_, 2*p_ + 1);
        ga_ij.set_val(k, wrt);
      }
      gradA_.set_val(i, j, ga_ij);
    }
  }
}


VecOfMats<cmplx>::type ASolver::calc_L()
{
  VecOfMats<cmplx>::type L (N_);\
  
  int i, j;
  MyMatrix<cmplx> inner, expand;
  for (i = 0; i < N_; i++)
  {
    inner = MyMatrix<cmplx> (p_, 2*p_ + 1);
    for (j = 0; j < N_; j++)
    {
      if (j == i) continue;
      expand = re_expandA(i, j);
      inner += expand;
    }
    L.set_val(i, inner);
  }
  return L;
}


void ASolver::print_Ei( int i, int p)
{
  cout << "This is my E for molecule " << i << " poles " << p <<  endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_E_ni(i,n,m).real();
      double im = get_E_ni( i, n, m).imag();
      r  = fabs( r) > 1e-9 ?  r : 0;
      im = fabs(im) > 1e-9 ? im : 0;
      cout << "(" << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}


void ASolver::print_Ai( int i, int p)
{
  cout << "This is my A for molecule " << i << " poles " << p <<  endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_A_ni(i,n,m).real();
      double im = get_A_ni( i, n, m).imag();
      r  = fabs( r) > 1e-9 ?  r : 0;
      im = fabs(im) > 1e-9 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}

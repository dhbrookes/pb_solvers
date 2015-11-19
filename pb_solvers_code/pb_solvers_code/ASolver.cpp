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
T_ (N_, N_), A_(N), reExpConsts_(sys.get_consts().get_kappa(),
                                 sys.get_lambda(), p), prevA_(N)
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

// one iteration of numerical solution for A
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
  
  x1 = expand_RX(  i, j);
  x2 = expand_SX(  i, j, x1);
  z  = expand_RHX( i, j, x2);
  return z;
  
//  int k, m;
//  cout << "This is my x1  for " << i << "  " << j << endl;
//  for(k = 0; k < p_; k++)
//  {
//    for(m = 0; m <= k; m++)
//    {
//      double  r = x1(k,m+p_).real();
//      double im = x1(k,m+p_).imag();
//      r  = fabs( r) > 1e-9 ?  r : 0;
//      im = fabs(im) > 1e-9 ? im : 0;
//      cout << " (" << r << "," << im << ")  ";
//    }
//    cout << endl;
//  }
//  
//  cout << "This is my x2  for " << i << "  " << j << endl;
//  for(k = 0; k < p_; k++)
//  {
//    for(m = 0; m <= k; m++)
//    {
//      double  r = x2(k,m+p_).real();
//      double im = x2(k,m+p_).imag();
//      r  = fabs( r) > 1e-9 ?  r : 0;
//      im = fabs(im) > 1e-9 ? im : 0;
//      cout << " (" << r << "," << im << ")  ";
//    }
//    cout << endl;
//  }
//  
//  cout << "This is my z for " << i << "  " << j << endl;
//  for(k = 0; k < p_; k++)
//  {
//    for(m = 0; m <= k; m++)
//    {
//      double  r = z(k,m+p_).real();
//      double im = z(k,m+p_).imag();
//      r  = fabs( r) > 1e-9 ?  r : 0;
//      im = fabs(im) > 1e-9 ? im : 0;
//      cout << " (" << r << "," << im << ")  ";
//    }
//    cout << endl;
//  }
}

// perform first part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RX(int i, int j)
{
  int n, m, s;
  cmplx inter;
  MyMatrix<cmplx> x1(p_, 2*p_ + 1);
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_(i,j).isSingular())
      {
        Pt vec = T_(i,j).get_TVec();
        if (vec.theta() > M_PI/2.0)
          inter = (n % 2 == 0 ? get_A_ni(j, n, -m) : -get_A_ni(j, n, -m));
        
        else
          inter = get_A_ni(j, n, m);
        
      } else
      {
        for (s = -n; s <= n; s++)
        {
          inter += T_(i, j).get_rval(n, m, s) * get_A_ni(j, n, s);
        } // end s
      }
      x1.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x1;
}

// perform second part of T*A and return results
MyMatrix<cmplx> ASolver::expand_SX(int i, int j, MyMatrix<cmplx> x1)
{
  int n, m, l;
  cmplx inter;
  MyMatrix<cmplx> x2(p_, 2*p_ + 1);
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = 0; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        inter += T_(i, j).get_sval(n, l, m) * x1(l, m+p_);
      } // end l
      
      x2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x2;
}

// perform third part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RHX(int i, int j, MyMatrix<cmplx> x2)
{
  int n, m, s;
  cmplx inter;
  MyMatrix<cmplx> z(p_, 2*p_ + 1);
  
  //fill zj:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (s = -n; s <= n; s++)
      {
        inter += conj(T_(i, j).get_rval(n, s, m)) * x2(n, s+p_);
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


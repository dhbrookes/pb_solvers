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
:p_(p), gamma_(N),delta_(N), E_(N), N_(sys.get_n()), T_ (N_, N_), A_(N),
reExpConsts_(sys.get_consts().get_kappa(), sys.get_lambda(), p), prevA_(N)
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
//  cout << "This is my R in solve_A" << endl;
//  for (int m = 0; m < 9; m++)
//  {
//    cout << "\t---m = " << m << "---" << endl;
//    for (int l = m; l < 9; l++)
//    {
//      for (int n = m; n <= l; n++)
//      {
//        double r = fabs(T_(0,1).get_rval(n, m, l).real())>1e-15 ?
//        T_(0,1).get_rval(n, m, l).real() : 0;
//        double im = fabs(T_(0,1).get_rval(n, m, l).imag())>1e-15 ?
//        T_(0,1).get_rval(n, m, l).imag() : 0;
//        cout << "(" << r << "," << im << ") | ";
//      }
//      cout << endl;
//    }
//  }
//  cout << endl;
  
//  cout << "This is my S in solve_A" << endl;
//  for (int m = 0; m < 9; m++)
//  {
//    cout << "\t---m = " << m << "---" << endl;
//    for (int l = m; l < 9; l++)
//    {
//      for (int n = m; n <= l; n++)
//      {
//        double r = fabs(T_(0,1).get_sval(n, l, m))>1e-15 ?
//                    T_(0,1).get_sval(n, l, m) : 0;
//        cout << r << " | ";
//      }
//      cout << endl;
//    }
//  }
//  cout << endl;
  
  prevA_ = A_;
  while(calc_change() > prec)
  {
    iter();
  }
  
  cout << endl;

}

// one iteration of numerical solution for A
void ASolver::iter()
{
  int i, j;
  MyMatrix<cmplx> Z, zj, ai;
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
  int i, k;
  double change = 0;
  cmplx prev, curr; // intermediate values
  for (i = 0; i < N_; i++)
  {
    for(k = 0; k < p_ * p_; k++)
    {
      prev = prevA_[i](k, k);
      curr = A_[i](k, k);
      change += abs(prev - curr) / abs(prev + curr);
    }
  }
  change *= 1 / (2*N_*p_*p_);
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

//re-expand element i of A with element (i, j) of T and return results
MyMatrix<cmplx> ASolver::re_expandA(int i, int j)
{
  
  int n, m, s, l;
  cmplx inter; // intermediate sum
  MyMatrix<cmplx> x1, x2, z;
      
  z = MyMatrix<cmplx>  (p_, 2*p_ + 1);
  x1 = MyMatrix<cmplx> (p_, 2*p_ + 1);
  x2 = MyMatrix<cmplx> (p_, 2*p_ + 1);
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (s = -n; s <= n; s++)
      {
        inter += T_(i, j).get_rval(n, m, s) * get_A_ni(i, n, s);
      } // end s
      x1.set_val(n, m+p_, inter);
    } // end m
  } //end n
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        inter += T_(i, j).get_sval(n, l, m) * x1(l, m+p_);
      } // end l
      x2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  
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
  d *= pow(ai, 2.0*n_dub+1.0) / (2.0*n_dub+1.0);
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
  for (i = 0; i < N_; i++)
  {
    A_.set_val(i, gamma_[i] * E_[i]);
  }
}


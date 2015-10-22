//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"
#include <iostream>

ASolver::ASolver(const int N, const int p, const BesselCalc* _bcalc,
                 SHCalc* _shCalc, System* sys)
:p_(p), _besselCalc_(_bcalc), _consts_(&sys->get_consts()), gamma_(N, N)
,delta_(N, N), E_(N), _shCalc_(_shCalc), _sys_(sys), N_(sys->get_n()),
T_ (N_, N_)
{
  // precompute all SH:
  all_sh.reserve(N_);
  int i = 0;
  int j = 0;

  // Calculate the T matrix (i.e the re expansion coefficients along
  // every inter-molecular vector (distance between every molecular center
  Pt v, ci, cj;  // inter molecular vector
  
//  all_inter_sh(N_, N_);
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      
      if (i == j) continue;
      ci = _sys_->get_centeri(i);
      cj = _sys_->get_centeri(j);
      v = ci - cj;
      // calculate spherical harmonics for inter molecular vector:
      _shCalc_->calc_sh(v.theta(), v.phi());
//      all_inter_sh.set_val(i, j, _shCalc_->get_full_result());
      T_.set_val(i, j, ReExpCoeffs(p_, v, _shCalc_->get_full_result(),
                                   _besselCalc_, _reExpConsts_,
                                   _consts_->get_kappa(),
                                   _sys_->get_lambda()));
    }
  }
  
  for (i = 0; i < N_; i++)
  {
    all_sh.push_back(calc_mol_sh(_sys_->get_molecule(i)));
  }
  
  //precompute gamma, delta and E:
  compute_gamma();
  compute_delta();
  compute_E();
  
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
  d *= pow(ai, 2.0*n_dub+1.0) / (2.0*n_dub+1.0);
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
  MyMatrix<double> gi;
  int i, j;
  for (i = 0; i < N_; i++)
  {
    gi = MyMatrix<double> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      gi.set_val(j, j, calc_indi_gamma(i, j));
    }
    gamma_.set_val(i, i, gi);
  }
}


/*
 Constructs the delta matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_delta()
 */
void ASolver::compute_delta()
{
  MyMatrix<double> di;
  int i, j;
  for (i = 0; i < N_; i++)
  {
    di = MyMatrix<double> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      di.set_val(j, j, calc_indi_delta(i, j));
    }
    delta_.set_val(i, i, di);
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
        ei.set_val(n, m + n, calc_indi_e(i, n, m));
      }
    }
    E_.set_val(i, ei);
  }
}


ASolver::~ASolver()
{
// potentially problematic. Need to address later
//  delete _besselCalc_;
//  delete _shCalc_;
}


ASolver::ASolver(const ASolver& other)
:_sys_(other._sys_), p_(other.p_), gamma_(other.gamma_),
delta_(other.delta_), N_(other.N_), E_(other.E_),
_consts_(other._consts_), all_sh(other.all_sh),
T_(other.T_), _besselCalc_(other._besselCalc_),
_reExpConsts_(other._reExpConsts_)
{
}

ASolver& ASolver::operator=(const ASolver& other)
{
  _besselCalc_ = other._besselCalc_;
  _shCalc_ = other._shCalc_;
  
  _sys_ = other._sys_;
  gamma_ = MatOfMats<double>::type(other.gamma_);
  delta_ = MatOfMats<double>::type(other.delta_);
  E_ = VecOfMats<cmplx>::type(other.E_);
  N_ = int(other.N_);
  p_ = int(other.p_);
  _consts_ =other._consts_;
  all_sh = vector<vector<MyMatrix<cmplx> > >(other.all_sh);
  return *this;
}


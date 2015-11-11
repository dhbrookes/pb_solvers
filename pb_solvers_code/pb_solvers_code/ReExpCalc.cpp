//
//  ReExpCalc.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 10/5/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ReExpCalc.h"

ReExpCoeffsConstants::ReExpCoeffsConstants(double kappa,
                                           double lambda, int p)
: kappa_(kappa), lambda_(lambda), p_(p),
a_(2*p, 4*p), b_(2*p, 4*p),
alpha_(2*p, 4*p), beta_(2*p, 4*p),
nu_(2*p, 4*p), mu_(2*p, 4*p)
{
  calc_a_and_b();
  calc_alpha_and_beta();
  calc_nu_and_mu();
}

void ReExpCoeffsConstants::calc_a_and_b()
{
  int m, n;
  double sign, a_val, b_val;
  
  //calculate a and b:
  for (n = 0; n < 2*p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      double nD = (double) n;
      double mD = (double) m;
      a_val = sqrt(((nD+mD+1.) * (nD-mD+1.)) / ((2.*nD+1.)* (2.*nD+3.)));
      if (m < 0)        sign = -1.0;
      else              sign = 1.0;
      b_val = sign * sqrt(((nD-mD-1.) * (nD-mD)) / ((2.*nD-1.) * (2.*nD+1.)));
      
      set_a_val(n, m, a_val);
      set_b_val(n, m, b_val);
    }
  }
}


void ReExpCoeffsConstants::calc_alpha_and_beta()
{
  int m, n;
  double alpha_val, beta_val;
  vector<double> alpha_m, beta_m;
  
  //calculate alpha and beta:
  for (n = 0; n < 2*p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      double nD = (double) n;
      double mD = (double) m;
      alpha_val = sqrt((nD + mD + 1.0) * (nD - mD + 1.0));
      beta_val  = (pow(lambda_*kappa_, 2.0)*alpha_val);
      beta_val *= (1.0 / ((2.0 * nD + 1.0) * (2.0 * nD + 3.0)));
      set_alpha( n, m, alpha_val);
      set_beta(  n, m, beta_val);
      
    }
  }
}


void ReExpCoeffsConstants::calc_nu_and_mu()
{
  int m, n, sign;
  double nu_val, mu_val;
  
  //calculate alpha and beta:
  for (n = 0; n < 2*p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      if (m < 0)        sign = -1.0;
      else              sign = 1.0;
      double nD = (double) n;
      double mD = (double) m;
      nu_val = sign * sqrt((nD - mD - 1.) * (nD - mD));
      mu_val = (pow(lambda_*kappa_, 2) * nu_val) / ((2.*nD-1.) * (2.*nD+1.));
      set_nu( n, m, nu_val);
      set_mu( n, m, mu_val);
    }
  }
}


ReExpCoeffs::ReExpCoeffs(Pt v, int p, MyMatrix<cmplx> Ytp,
                               BesselCalc bessel_calc,
                               ReExpCoeffsConstants _consts,
                               double kappa, double lambda)
: v_(v), p_(p), Ytp_(Ytp), _besselCalc_(bessel_calc),
kappa_(kappa), lambda_(lambda), _consts_(_consts)
{
  
  if (_besselCalc_.get_num_vals() < 2 * p_)
  {
    throw BesselSizeException(p_, _besselCalc_.get_num_vals());
  }
  
  if (Ytp_.get_nrows() < 2 * p_)
  {
    throw SHSizeException(p_, Ytp_.get_nrows());
  }
  calc_r();
  calc_s();
}


void ReExpCoeffs::calc_r()
{
  int n, m, s;
  cmplx val;
  cmplx ic = cmplx(0, 1);
  double phi = v_.phi();
  double theta = v_.theta();
  R_ = MyVector<MyMatrix<cmplx> > (2*p_); // n range of 0 to 2p-1 is needed!
  for (n = 0; n < 2 * p_; n++)
  {
    R_.set_val(n, MyMatrix<cmplx> (2*p_, 4*p_));//s range: -2p+1 to 2p-1 needed!
    for (s = -n; s <= n; s++)
    {
      val = get_yval(n, -s);
      set_rval(n, 0, s, val);
    }
  }
  
  for (m = 0; m < p_; m++)
  {
    for (n=m+2; n < 2 * p_ - m; n++)
    {
      for (s=-n+1; s < n; s++)
      {
        val = 0.5 * exp(-ic * phi) * (1 + cos(theta)) *
        _consts_.get_b_val(n, s-1) * get_rval(n, m, s-1);
        val -= 0.5 * exp(ic * phi) * (1 - cos(theta)) *
        _consts_.get_b_val(n, -s-1) * get_rval(n, m, s+1);
        val += sin(theta) * _consts_.get_a_val(n ,s) *
        get_rval(n, m, s);
        val *= ( 1.0 / _consts_.get_b_val( n, m) );
        set_rval(n-1, m+1, s, val);
      }
    }
  }
}

void ReExpCoeffs::calc_s()
{
  int m, n, l;
  double val;
  double r = v_.r();
  S_ = MyVector<MyMatrix<double> > ( 2 * p_ );
  vector<double> besselK = _besselCalc_.calc_mbfK( 2*p_, kappa_*r);
  
  for (n = 0; n < 2*p_; n++)
  {
    S_.set_val(n, MyMatrix<double> ( 2*p_, 4*p_));
  }
  
  for (l = 0; l < 2 * p_; l++)
  {
    val = pow((lambda_ / r), l) * (besselK[l] * exp(-kappa_*r)) / r;
    set_sval( 0, l, 0, val );
    
    // Set S_{n,0}^0 vals
    set_sval( l, 0, 0, pow(-1.0, l) * val );
  }
  
  for (l = 1; l < 2 * p_ - 2; l++)
  {
    val  = _consts_.get_beta(l-1, 0) * get_sval( 0, l-1, 0);
    val += _consts_.get_alpha( l, 0) * get_sval( 0, l+1, 0);
    val *= -1.0 / _consts_.get_alpha(0, 0);
    set_sval( 1, l, 0, val );
  }
  
  for (n = 1; n < p_ - 1; n++)
  {
    for(l = n + 1; l < 2*p_ - n - 2; l++)
    {
      val  = _consts_.get_beta(l-1, 0) * get_sval(  n, l-1, 0);
      val += _consts_.get_beta(n-1, 0) * get_sval(n-1,   l, 0);
      val += _consts_.get_alpha( l, 0) * get_sval(  n, l+1, 0);
      val *= -1.0 / _consts_.get_alpha(n, 0);
      set_sval(n+1, l, 0, val);
    }
  }
  
  /* old code implementation differs from paper:
   for (int m = 1; m < m_p; m++)
   {
   for (int l = m; l < 2*m_p-1-m; l++)	// EQ 1.7, Lotan 2006
   U[l][m][m] = -(DELTA(l,-m)*U[l-1][m-1][m-1] +
   GAMMA(l+1,m-1)*U[l+1][m-1][m-1])/GAMMA(m,-m);
   
   for (int n = m; n < m_p-1; n++)				// EQ 1.8, Lotan 2006
   for (int l = n+1; l < 2*m_p-2-n; l++)
   U[l][n+1][m] = -(BETA(l-1,m)*U[l-1][n][m] +
   BETA(n-1,m)*U[l][n-1][m] +
   ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
   }
   */
  double val2;
  for (m = 1; m < p_ ; m++)
  {
    for (l = m; l < 2*p_ - m - 1; l++)
    {
      val  = _consts_.get_mu(  l,  -m) * get_sval(m-1, l-1, m-1);
      val += _consts_.get_nu(l+1,m-1) * get_sval( m-1, l+1, m-1);
      val *= ( -1.0 / _consts_.get_nu( m, -m) );
      set_sval(m, l, m, val);
    }
    
    for (n = m; n < p_ - 1; n++)
    {
      for (l = n + 1; l < 2*p_ - n - 2; l++)
      {
        val2  = _consts_.get_beta(l-1, m) * get_sval(  n, l-1, m);
        val2 += _consts_.get_beta(n-1, m) * get_sval(n-1,   l, m);
        val2 += _consts_.get_alpha( l, m) * get_sval(  n, l+1, m);
        val2 *= (-1.0 / _consts_.get_alpha(n, m));
        set_sval(n+1, l, m, val2);
      }
    }
  }
  
  // Filling in the rest !
  for (n = 1; n < p_ ; n++)
    for (l = 0; l < p_; l++)
      for (m = -n; m <= n; m++)
        if ( n  > l )
          set_sval(n, l, m, pow(-1.0, n) * get_sval( l, n, m));
  
  
  for (n = 1; n < p_ ; n++)
    for (l = 0; l < p_; l++)
      for (m = -n; m <= n; m++)
        if ( m  < 0 )
          set_sval(n, l, m, get_sval( l, n, -m));
  
  
} // end calc_s


ReExpCoeffs& ReExpCoeffs::operator=(ReExpCoeffs& other)
{
  R_ = other.R_;
  S_ = other.S_;
  
  v_ = other.v_;
  p_ = other.p_;
  Ytp_ = other.Ytp_;
  _besselCalc_ = other._besselCalc_;
  _consts_ = other._consts_;
  
  kappa_ = other.kappa_;
  lambda_ = other.lambda_;
  
  return *this;
}




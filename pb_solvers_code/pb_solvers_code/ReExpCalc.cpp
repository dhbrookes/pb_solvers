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
:p_(p), a_(4*p, 4*p), b_(4*p, 4*p)
{
  calc_a_and_b();
 // calc_alpha_and_beta();
 // calc_nu_and_mu();
}

void ReExpCoeffsConstants::calc_a_and_b()
{
  
  int m, n,sign;
  double a_val, b_val;
  
  //calculate a and b:
  for (n = 0; n <= 2*p_; n++)
  {
    for (m = -n+1; m < n; m++)
    {
      if (n < (m-2))
      {
        a_val = 0.0;
        b_val = 0.0;
      }
      else
      {
        a_val = sqrt(((n+m+1) * (n-m+1)) / ((2*n+1)* (2*n+3)));
        if (m < 0)        sign = -1.0;
        else if (m == 0)  sign = 0.0;
        else              sign = 1.0;
        b_val = sign * sqrt(((n-m-1) * (n-m)) / ((2*n-1) * (2*n+1)));
      }
      set_a_val(m, n, a_val);
      set_b_val(m, n, b_val);
    }
  }
}
//
//
//void ReExpCoeffsConstants::calc_alpha_and_beta()
//{
//
//
//  int m, n;
//  double alpha_val, beta_val;
//  vector<double> alpha_m, beta_m;
//  
//  //calculate alpha and beta:
//  for (m = 0; m < p_; m++)
//  {
//    for (n = -1; n < 2*p_; n++)
//    {
//      alpha_val = sqrt((n + m + 1) * (n - m + 1));
//      beta_val = (pow(lambda_, 2)*kappa_*alpha_val) / ((2*n+1)*(2*n+3));
//      alpha_.set_val(m, n, alpha_val);
//      beta_.set_val(m, n, beta_val);
//    }
//  }
//}
//
//
//void ReExpCoeffsConstants::calc_nu_and_mu()
//{
//  
//  int m, n, sign;
//  double nu_val, mu_val;
//  
//
//  //calculate alpha and beta:
//  for (m = -p_+1; m < p_; m++)
//  {
//    for (n = -1; n < 2*p_; n++)
//    {
//      if (m < 0)        sign = -1.0;
//      else if (m == 0)  sign = 0.0;
//      else              sign = 1.0;
//      nu_val = sign * sqrt((n - m - 1) * (n - m));
//      mu_val = (pow(lambda_, 2) * kappa_ * nu_val) / ((2*n-1) * (2*n+1));
//      
//    }
//  }
//}


ReExpCoeffs_IJ::ReExpCoeffs_IJ(int p, ShPt v, MyMatrix<cmplx>* Ytp,
                               BesselCalc * BesselCalc,
                               ReExpCoeffsConstants* _consts,
                               double kappa, double lambda)
:p_(p), v_(v), _Ytp_(Ytp), _besselCalc_(BesselCalc),
 kappa_(kappa), lambda_(lambda), _consts_(_consts)
{
  if (_besselCalc_->get_num_vals() < 2 * p_)
  {
   throw BesselSizeException(p_, _besselCalc_->get_num_vals());
  }
  
  if (_Ytp_->get_nrows() < 2 * p_)
  {
    throw SHSizeException(p_, _Ytp_->get_nrows());
  }
  calc_r();
  calc_s();
}


void ReExpCoeffs_IJ::calc_r()
{
  int n, m, s;
  cmplx val;
  cmplx ic = cmplx(0, 1);
  double phi = v_.get_phi();
  double theta = v_.get_theta();
  R_ = MyVector<MyMatrix<cmplx> > (2*p_); // n range of 0 to 2p-1 is needed!
  for (n = 0; n < 2 * p_; n++)
  {
    R_.set_val(n, MyMatrix<cmplx> (2*p_, 4*p_)); // s range: -2p+1 to 2p-1 needed!
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
          _consts_->get_b_val(s-1, n) * get_rval(n, m, s-1);
        val -= 0.5 * exp(ic * phi) * (1 - cos(theta)) *
          _consts_->get_b_val(-s- 1, n) * get_rval(n, m, s+1);
        val += sin(theta) * _consts_->get_a_val(s, n) *
          get_rval(n, m, s);
        val *= 1 / _consts_->get_b_val(m, n);
        set_rval(n-1, m+1, s, val);
      }
    }
  }
}

void ReExpCoeffs_IJ::calc_s()
{
  int m, n, l;
  cmplx val;
  double r = v_.get_r();
  S_ = MyVector<MyMatrix<cmplx> > ( 2 * p_ );
  vector<double> besselK = _besselCalc_->calc_mbfK( 2*p_, kappa_*r);
  
  for (n = 0; n < 2*p_; n++)
  {
    S_.set_val(n, MyMatrix<cmplx> ( 2*p_, 4*p_));
  }
  
  for (l = 0; l < 2 * p_; l++)
  {
    val = pow((lambda_ / r), l) * (besselK[l] * exp(-kappa_*r)) / r;
    set_sval( 0, l, 0, val );
    
    // Set S_{n,0}^0 vals
    set_sval( l, 0, 0, pow(-1.0, l) * val );
    }
  
  for (n = 1; n < 2*p_ - 1; n++)
  {
    for(l = n+1; l < 2*p_ - n - 1; l++)
    {
      val  = calc_beta(0, l-1) * get_sval( n, l-1, 0);
      val += calc_beta(0, n-1) * get_sval( n,   l, 0);
      val += calc_alpha(0,  l) * get_sval( n, l+1, 0);
      val *= -1.0 / calc_alpha(0, n);
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
  cmplx val2;
  for (m = 1; m < p_ - 1; m++)
  {
    for (l = m; l < 2*p_ - m - 1; l++)
    {
      val  = calc_mu( -m,   l) * get_sval(m-1, l-1,m-1);
      val += calc_nu(m-1, l+1) * get_sval(m-1, l+1,m-1);
      val *= 1.0 / calc_nu(-m, m);
      set_sval(m+1, l, m+1, val);
    }
    
    for (n = m; n < p_ - 1; n++)
    {
      for (l = n + 1; l < 2*p_ - n - 1; l++)
      {
        val2  = calc_beta(m, l-1) * get_sval(  n, l-1, m);
        val2 += calc_beta(m, n-1) * get_sval(n-1,   l, m);
        val2 += calc_alpha(m,  l) * get_sval(  n, l+1, m);
        val2 *= -1.0 / calc_alpha(m, n);
        set_sval(n+1, l, m+1, val2);
      }
    }
  }
  
} // end calc_s


//const double ReExpCoeffs_IJ::calc_a(int m, int n)
//{
//  double mD = (double) m;
//  double nD = (double) n;
//  double a_val = sqrt(((nD+mD+1.0) * (nD-mD+1.0)) /
//                      ((2.0*nD+1.0)* (2.0*nD+3.0)));
//  return a_val;
//}
//
///*
// 
// Original code: 0 is positive
// for (int n = 0; n < 2*N_POLES; n++)
//   for (int m = 0; m <= n; m++)
//   {
//     b(n,m) = sqrt((REAL)(n-m-1)*(n-m)/((2*n+1)*(2*n-1)));//Lotan 2006 (eq 1.2)
//     if (m != 0)
//       b(n,-m) = -sqrt((REAL)(n+m-1)*(n+m)/((2*n+1)*(2*n-1)));
//   }
// */
//const double ReExpCoeffs_IJ::calc_b(int m, int n)
//{
//  double sign = 1.0;
//  if (m < 0)  sign = -1.0;
//  double mD = (double) m;
//  double nD = (double) n;
// 
//  double b_val = sign * sqrt(((nD-mD-1.0) * (nD-mD)) /
//                             ((2.0*nD-1.0) * (2.0*nD+1.0)));
//  return b_val;
//}


const double ReExpCoeffs_IJ::calc_alpha(int m, int n)
{
  double mD = (double) m;
  double nD = (double) n;
  double alpha_val = sqrt((nD + mD + 1.0) * (nD - mD + 1.0));
  return alpha_val;
}


const double ReExpCoeffs_IJ::calc_beta(int m, int n)
{
  double alpha_val = calc_alpha(m, n);
  double beta_val = (pow(lambda_, 2)*kappa_*alpha_val) / ((2*n+1)*(2*n+3));
  return beta_val;
}


const double ReExpCoeffs_IJ::calc_nu(int m, int n)
{
  int sign;
  if (m < 0)        sign = -1.0;
  else if (m == 0)  sign = 0.0;
  else              sign = 1.0;
  double nu_val = sign * sqrt((n - m - 1) * (n - m));
  return nu_val;
}


const double ReExpCoeffs_IJ::calc_mu(int m, int n)
{
  double mu_val = (pow(lambda_, 2) * kappa_ * calc_nu(m, n))
                  / ((2*n-1) * (2*n+1));
  return mu_val;
}






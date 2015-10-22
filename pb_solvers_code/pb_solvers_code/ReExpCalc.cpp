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
:p_(p), a_(p, p), b_(p, p), alpha_(p, p+1), beta_(p, p+1)
,nu_(2*p, p), mu_(2*p, p)
{
  calc_a_and_b();
  calc_alpha_and_beta();
  calc_nu_and_mu();
}

void ReExpCoeffsConstants::calc_a_and_b()
{
  
  int m, n,sign;
  double a_val, b_val;
  
  //calculate a and b:
  for (m = 0; m < p_; m++)
  {
    for (n = 0; n < 2*p_-m; n++)
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
      a_.set_val(m, n, a_val);
      b_.set_val(m, n, b_val);
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


ReExpCoeffs::ReExpCoeffs(int p, Pt v, MyMatrix<cmplx>* Ytp,
                               const BesselCalc* BesselCalc,
                               ReExpCoeffsConstants* _consts,
                               double kappa, double lambda)
:p_(p), v_(v), Ytp_(Ytp), kappa_(kappa), lambda_(lambda)
{
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
    R_.set_val(n, MyMatrix<cmplx> (p_, 2*p_));
    for (s = -n; s <= n; n++)
    {
      val = Ytp_->operator()(n, -s);
      R_[n].set_val(0, s, val);
      
    }
  }
  
  for (m = 0; m < p_; m++)
  {
    for (n=m+2; n< 2*p_ - m; n++)
    {
      for (s=-n+1; s < n; s++)
      {
        val = 0.5 * exp(ic * phi) * (1 + cos(theta)) *
        calc_b(s-1, n) * R_[n](m, s-1);
        val -= 0.5 * exp(ic * phi) * (1 - cos(theta)) *
        calc_b(-s+1, n) * R_[n](m, s-1);
        val += sin(theta)*calc_a(s, n)*R_[n](m, s);
        val *= 1 / calc_b(m, n);
        R_[n-1].set_val(m+1, s, val);
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
  vector<double> besselK = _besselCalc_->calc_mbfK( 2*p_, kappa_*r);
  
  for (n = 0; n < 2*p_; n++)
  {
    S_.set_val(m, MyMatrix<cmplx> (2*p_, 2*p_));
  }
  
  for (l = 0; l < 2*p_; l++)
  {
    val = pow((lambda_ / r), l) * (besselK[l-1] * exp(-kappa_*r)) / r;
    S_[0].set_val(0, l, val);
  }
  
  for (n = 0; p_ - 1; n++)
  {
    for(l = n+1; l < 2*p_ - n - 1; l++)
    {
      val = calc_beta(0, l-1) * S_[0](n, l-1);
      val += calc_beta(0, n-1) * S_[0](n-1, l);
      val += calc_alpha(0, l) * S_[0](n, l+1);
      val *= -1 / calc_alpha(0, n);
      S_[0].set_val(n+1, l, val);
    }
  }
  
  cmplx val2;
  int l2;
  for (m=0; m < p_-1; m++)
  {
    for (l = m; l < 2*p_ - m - 1; l++)
    {
      val = calc_mu(-m-1, l) * S_[m](m, l-1);
      val += calc_nu(m, l+1) * S_[m](m, l+1);
      val *= 1 / calc_nu(-m-1, m+1);
      S_[m+1].set_val(m+1, l, val);
      
      for (n = m; n < p_-1; n++)
      {
        for (l2 = n+1; l2 < 2*p_ - n -1; l2++)
        {
          val2 = calc_beta(m+1, l2-1) * S_[m+1](n, l-1);
          val2 += calc_beta(m+1, n-1) * S_[m+1](n-1, l);
          val2 += calc_alpha(m+1, l) * S_[m+1](n-1, l);
          val2 *= -1/(calc_alpha(m+1, n));
          S_[m+1].set_val(n+1, l, val2);
        }
      }
      
    }
  }
  
}

const double ReExpCoeffs_IJ::calc_a(int m, int n)
{
  double a_val = sqrt(((n+m+1) * (n-m+1)) / ((2*n+1)* (2*n+3)));
  return a_val;
}



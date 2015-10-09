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
:lambda_(lambda), p_(p), kappa_(kappa)
{
  calc_a_and_b();
  calc_alpha_and_beta();
  calc_nu_and_mu();
}

void ReExpCoeffsConstants::calc_a_and_b()
{
  a_.reserve(p_-1);
  b_.reserve(p_-1);
  
  int m, n, inner_size, sign;
  double a_val, b_val;
  vector<double> a_m, b_m;
  
  //calculate a and b:
  for (m = 0; m < p_-1; m++)
  {
    inner_size = (2*p_-m-1);
    a_m.reserve(inner_size);
    b_m.reserve(inner_size);
    for (n = 0; n < 2*p_-m-1; n++)
    {
      if (n < (m-2))
      {
        a_val = 0.0;
        b_val = 0.0;
      }
      else
      {
        a_val = sqrt(((n+m+1) * (n-m+1)) / ((2*n+1) * (2*n+3)));
        if (m < 0)        sign = -1.0;
        else if (m == 0)  sign = 0.0;
        else              sign = 1.0;
        b_val = sign * sqrt(((n-m-1) * (n-m)) / ((2*n-1) * (2*n+1)));
      }
      a_m.push_back(a_val);
      b_m.push_back(b_val);
    }
    a_.push_back(a_m);
    b_.push_back(b_m);
    a_m.clear();
    b_m.clear();
  }
}


void ReExpCoeffsConstants::calc_alpha_and_beta()
{
  alpha_.reserve(p_-1);
  beta_.reserve(p_-1);

  int m, n, inner_size;
  double alpha_val, beta_val;
  vector<double> alpha_m, beta_m;
  
  
  //calculate alpha and beta:
  for (m = 0; m < p_-1; m++)
  {
    inner_size = 2*p_-1;
    alpha_m.reserve(inner_size);
    beta_m.reserve(inner_size + 1);
    for (n = -1; n < 2*p_-1; n++)
    {
      alpha_val = sqrt((n + m + 1) * (n - m + 1));
      beta_val = (pow(lambda_, 2)*kappa_*alpha_val) / ((2*n + 1)*(2*n+3));
      alpha_m.push_back(alpha_val);
      beta_m.push_back(beta_val);
    }
    alpha_.push_back(alpha_m);
    beta_.push_back(beta_m);
    alpha_m.clear();
    beta_m.clear();
  }
}


void ReExpCoeffsConstants::calc_nu_and_mu()
{
  nu_.reserve(2*(p_-1));
  mu_.reserve(2*(p_-1));
  
  int m, n, inner_size, sign;
  double nu_val, mu_val;
  vector<double> nu_m, mu_m;
  

  //calculate alpha and beta:
  for (m = -p_+1; m < p_-1; m++)
  {
    inner_size = 2*p_-1;
    nu_m.reserve(inner_size);
    mu_m.reserve(inner_size);
    for (n = -1; n < 2*p_-1; n++)
    {
      if (m < 0)        sign = -1.0;
      else if (m == 0)  sign = 0.0;
      else              sign = 1.0;
      nu_val = sign * sqrt((n - m - 1) * (n - m));
      mu_val = (pow(lambda_, 2) * kappa_ * nu_val) / ((2*n-1) * (2*n+1));
      nu_m.push_back(nu_val);
      mu_m.push_back(mu_val);
    }
    nu_.push_back(nu_m);
    mu_.push_back(mu_m);
    nu_m.clear();
    mu_m.clear();
  }
}



//
//  BesselCalc.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "BesselCalc.h"


BesselConstants::BesselConstants(const int N)
:numVals_(N)
{
  kConsts_.reserve(numVals_);
  int n;
  double kval;
  for (n = 0; n < numVals_; n++)
  {
    kval = 1.0 / ((2*n+1) * (2*n-1));
    kConsts_.push_back(kval);
  }
}


BesselCalc::BesselCalc(int N, BesselConstants* _consts)
: numVals_(N)
{
  assert (_consts->get_n() == numVals_);
  _consts_ = new BesselConstants;
  _consts_ = _consts;
}

/*
 Modified Bessel K function calculated using the recursion:
 k n+1(z) = k n(z) + ( k n-1(z) * z^2 )/( (2n+1)*(2n-1) )
 
 param n is the number of bessel functions to calculate
  Will return vector for \hat{i} that ranges from 0 to nval-1
 */
const vector<double> BesselCalc::calc_mbfK(const int n,
                                           const double z) const
{
  vector<double> K;
  K.reserve(n);
  
  if (n > 0) K.push_back(1.0);
  if (n > 1) K.push_back(1 + z);
  
  double z_sq = z * z;
  int i;
  double val;
  for (i = 2; i < n; i++)
  {
    val = K[i-1] + (z_sq * K[i-2] * _consts_->get_kconst_val(i-1));
    K.push_back(val);
  }
  return K;
}

/*
 Modified Bessel I function calculated using the recursion:
 i n(z) = 1 + sum_(j=1)^L t_j^n(z^2/2)
 where: t_j^n(y) = (1/j)*t_(j-1)^n(y)*(y/(2n+2j+3))
 
 param n is the number of bessel functions to calculate
     Will return vector for \hat{i} that ranges from 0 to nval-1
 */
const vector<double> BesselCalc::calc_mbfI(const int n,
                     const double z) const
{
  vector<double> I;
  I.reserve(n);
  for (int k = 0; k < n; k++)
  {
    I.push_back(1);
  }
  
  if (z != 0)
  {
    double y = 0.5 * z * z;
    int k, j;
    double t;
    for (k = 0; k < n; k++)
    {
      t = y / (2*k + 3);
      for (j = 1; j <= 20; j++)
      {
          I[k] += t;
        t *= y / ((j+1) * (2 * (k+j) + 3 ));  //EQ 1.15
        if (t < 1e-20) break;
      }
    }
  }
  return I;
}

BesselCalc::BesselCalc()
:numVals_(Constants::MAX_NUM_POLES)
{
  _consts_ = new BesselConstants;
  
}

BesselCalc::~BesselCalc()
{
//  delete _consts_;
}

BesselCalc::BesselCalc(const BesselCalc& other)
:numVals_(other.numVals_)
{
  _consts_ = new BesselConstants;
  _consts_ = other._consts_;
}

BesselCalc& BesselCalc::operator=(const BesselCalc& other)
{
  numVals_ = int(other.numVals_);
  _consts_ = new BesselConstants;
  _consts_ = other._consts_;
  return *this;
}





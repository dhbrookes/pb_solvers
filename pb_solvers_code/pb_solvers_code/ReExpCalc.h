//
//  ReExpCalc.h
//  pb_solvers_code
//
//  Created by David Brookes on 10/5/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef ReExpCalc_h
#define ReExpCalc_h

#include <stdio.h>
#include "MyMatrix.h"
#include "Constants.h"
#include "util.h"

using namespace std;

/*
 Class for pre-computing the constants in the re-expansion coefficne
 */
class ReExpCoeffsConstants
{
protected:
  vector<vector<double> > a_;  // first index is m, second is n
  vector<vector<double> > b_;
  vector<vector<double> > alpha_;
  vector<vector<double> > beta_;
  vector<vector<double> > nu_;
  vector<vector<double> > mu_;
  double lambda_; // uniform scaling factor (section 4.5 of Lotan 2006)
  int p_;
  double kappa_; //from Constants
  
  void calc_a_and_b();
  void calc_alpha_and_beta();
  void calc_nu_and_mu();

public:
  
  ReExpCoeffsConstants(double kappa, double lambda,
                       int p=Constants::MAX_NUM_POLES);
  
  const double get_a_val(int m, int n) const      { return a_[m][n]; }
  const double get_b_val(int m, int n) const      { return b_[m][n]; }
  const double get_alpha_val(int m, int n) const  { return alpha_[m][n+1]; }
  const double get_beta_val(int m, int n) const   { return beta_[m][n+1]; }
  const double get_nu_val(int m, int n) const     { return nu_[m+(p_-1)][n]; }
  const double get_mu_val(int m, int n) const     { return mu_[m+(p_-1)][n];}
  
};

/*
 Class containing necessary information for calculating re-expansion
 coefficients (labeled T matrix in Lotan 2006)
 */
class ReExpCoeffs
{
protected:
  
  int p_; // max value of n when solving for A
  vector<MyMatrix<cmplx> >  R_;  // rotation coefficients
  vector<MyMatrix<cmplx> >  S_;  // translation coefficients
  ReExpCoeffsConstants*     _theseConsts_;
  Constants*                _consts_;
  
  void calc_r(double theta);  // calculate all the values for R_
  void calc_s(); // calculate all the values for S_
  
public:
  
  ReExpCoeffs();
  ReExpCoeffs(int p, ReExpCoeffsConstants* _consts);
  virtual ~ReExpCoeffs();
  ReExpCoeffs& operator=(const ReExpCoeffs* other);
  
  cmplx get_rval(int n, int m, int s) { return R_[n](m, s); }
  cmplx get_sval(int i, int n, int m) { return S_[i](n, m); }
  
  
};

#endif /* ReExpCalc_hpp */

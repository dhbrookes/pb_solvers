
//  ReExpCalc.h
//  pb_solvers_code
//
//  Created by David Brookes on 10/5/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef ReExpCalc_h
#define ReExpCalc_h

#include "Constants.h"
#include "BesselCalc.h"
#include "SHCalc.h"
#include "MyMatrix.h"
#include "util.h"

#include <sstream>


using namespace std;

class BesselSizeException
{
protected:
int p_;
int besselSize_;

public:
BesselSizeException(const int p, const int besselSize)
:p_(p), besselSize_(besselSize)
{
}

virtual const char* what() const throw()
{
ostringstream ss;
ss << "The bessel vector is the wrong size. It is supposed to be: " <<
2 * p_ <<" but is:  " << besselSize_ << endl;
return ss.str().c_str();
}

};


class SHSizeException
{
protected:
  int p_;
  int shSize_;
  
public:
  SHSizeException(const int p, const int shSize)
  :p_(p), shSize_(shSize)
  {
  }
  
  virtual const char* what() const throw()
  {
    ostringstream ss;
    ss << "The spherical harmonics vector is the wrong size." <<
    "It is supposed to be: " << 2 * p_ <<" but is:  " << shSize_ << endl;
    return ss.str().c_str();
  }
  
};


/*
 Class for pre-computing the constants in the re-expansion coefficne
 */
class ReExpCoeffsConstants
{
protected:
  MyMatrix<double> a_;  // first index is m, second is n
  MyMatrix<double> b_;
  MyMatrix<double> alpha_;
  MyMatrix<double> beta_;
  MyMatrix<double> nu_;
  MyMatrix<double> mu_;
  double lambda_; // uniform scaling factor (section 4.5 of Lotan 2006)
  int p_;
  double kappa_; //from Constants
  
  void calc_a_and_b();
  void calc_alpha_and_beta();
  void calc_nu_and_mu();
  
public:
  
  /*
   Initialize proper amount of memory in default constructor:
   */
  ReExpCoeffsConstants(double kappa, double lambda,
                       int p=Constants::MAX_NUM_POLES);
  
  double get_a_val(int n, int m)  { return a_(n, m + 2*p_); }
  double get_b_val(int n, int m)  { return b_(n, m + 2*p_); }
  double get_alpha(int n, int m)  { return alpha_(n, m + 2*p_); }
  double get_beta(int n, int m)   { return  beta_(n, m + 2*p_); }
  double get_nu(int n, int m)     { return nu_(n, m + 2*p_);}
  double get_mu(int n, int m)     { return mu_(n, m + 2*p_);}
  
  void set_a_val( int n, int m, double val) {a_.set_val(n, m + 2*p_, val);}
  void set_b_val(int n, int m, double val)  {b_.set_val(n, m + 2*p_, val);}
  void set_alpha( int n, int m, double val) {alpha_.set_val(n, m + 2*p_, val);}
  void set_beta(int n, int m, double val)   {beta_.set_val(n, m + 2*p_, val);}
  void set_nu( int n, int m, double val)    {nu_.set_val(n, m + 2*p_, val);}
  void set_mu(int n, int m, double val)     {mu_.set_val(n, m + 2*p_, val);}
  
};

/*
 Class representing one entry in the re-expansion coefficient matrix. So if
 that matrix is T (as in Lotan 2006), then this class contains the info
 for one T^(i,j)
 */
class ReExpCoeffs
{
protected:
  
  int p_; // max value of n when solving for A
  
  /*
   R_ contains rotation coefficients for this entry. R_ has three
   indices: R[n](m, s)
   And the range of each:  0 <= n <  poles
   -n <= m <= n, but -m = conj(+m), so really just [0,p)
   -n <= s <= n
   */
  MyVector<MyMatrix<cmplx> >  R_;
  /*
   S_ contains translation coefficients for this entry. S_ has three
   indices: S[m](n, l)
   And the range of each:  0 <= n <  poles
   0 <= l <= poles
   -n <= m <= n
   */
  MyVector<MyMatrix<double> >  S_;
  //  ReExpCoeffsConstants*     _consts_;
  double kappa_; //from Constants
  double lambda_; // uniform scaling factor (section 4.5 of Lotan 2006)
  
  Pt v_; //computing re-expansion along this vector
  ReExpCoeffsConstants* _consts_;
  
  /*
   Spherical harmonics for this v_:
   */
  MyMatrix<cmplx> Ytp_;
  const BesselCalc* _besselCalc_;
  
  void calc_r();  // calculate all the values for R_
  void calc_s(); // calculate all the values for S_
  
public:
  
  ReExpCoeffs() { };
  ReExpCoeffs(int p, Pt v, MyMatrix<cmplx> Ytp, const BesselCalc * BesselCalc,
                 ReExpCoeffsConstants* _consts, double kappa, double lambda);
  //  virtual ~ReExpCoeffs_IJ();
  //  ReExpCoeffs_IJ& operator=(const ReExpCoeffs_IJ* other);
  
  cmplx get_yval(int n, int s)
  {
    if ( s < 0 ) return conj(Ytp_(n, -s));
    else         return Ytp_(n, s);
  }
  
  cmplx get_rval(int n, int m, int s)
  {
    if ( m < 0 ) return conj(R_[n](-m, -s+2*p_));
    else         return R_[n](m, s+2*p_);
  }
  
  double get_sval(int n, int l, int m)
  { return S_[n](l, m+2*p_); }
  
  void set_rval(int n, int m, int s, cmplx val)
  {
    if ( m < 0 ) (&R_[n])->set_val(-m, -s+2*p_, conj( val ));
    else         (&R_[n])->set_val( m, s+2*p_, val );
  }
  
  void set_sval(int n, int l, int m, double val)
  { (&S_[n])->set_val(l, m+2*p_, val); }
  
};

#endif /* ReExpCalc_hpp */
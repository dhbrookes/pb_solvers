//
//  SHCalc.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "SHCalc.h"


SHCalcConstants::SHCalcConstants(const int N)
:nPoles_(N), legConsts1_(2*N, 2*N), legConsts2_(2*N, 2*N),
shConsts_(2*N, 2*N), dubFac_(2*N)
{
    vector<double> temp;
    temp.reserve(4 * nPoles_);
    temp.push_back(0);
    int i, n, m;
    for (i = 1; i < 4 * nPoles_; i++)
    {
        temp.push_back(temp[i-1] * sqrt(i));
    }

    for (n = 0; n < 2 * nPoles_; n++)
    {
        for (m = 0; m <= n; m++)
        {
            legConsts1_.set_val(n, m, (2*n-1) / (n-m));
            legConsts2_.set_val(n, m, (n+m-1) / (n-m));
            shConsts_.set_val(n, m, temp[n-m] / temp[n+m]);
        }
    }

    dubFac_[0] = 1.0;
    dubFac_[1] = 1.0;
    for (i = 2; i < 2 * nPoles_; i++)
    {
        dubFac_[i] = dubFac_[i-1] * (2*i - 1);
    }
    
}

/*
 Full calculation is performed in constructor:
 */
SHCalc::SHCalc(const int N, const SHCalcConstants* _consts,
               const double theta, const double phi)
:_consts_(_consts), nPoles_(N), P_(2 * nPoles_, 2 * nPoles_), Y_(2 * nPoles_, 2 * nPoles_),
theta_(theta), phi_(phi)
{
    assert (_consts_->get_n() == nPoles_);
    
    calc_legendre();
    calc_sh();
}

/*
 
Calculate the Legendre polynomial for the input theta using the
 recursion functions for the polynomials, which are as follows:

Pl,l (x) = (-1)^l * (2l-1)!! * (1-x^2)^(l/2)                            (1)
Pl,l+1 (x) = x * (2l+1) * Pl,l(x)                                       (2)
Pl,m (x) = x * (2l-1)/(l-m) * Pl-1,m(x) - (l+m-1)/(l-m) * Pl-2,m(x)     (3)
*/
void SHCalc::calc_legendre()
{
    double x = cos(theta_);
    P_.set_val(0, 0, 1.0);  // base value for recursion
    P_.set_val(1, 0, x);
    
    int l, m, lInd, mInd;
    double val;
    for (l = 0; l < 2 * nPoles_; l++)
    {
        for (m = 0; m < l; m++)
        {
            if ((l == 0 && m == 0) || (l == 1 && m == 0)) continue;  //skip the base values
            else if (l == m)
            {
                val = pow(-1, (double) l) * _consts_->get_dub_fac_val(l) * pow(1-x*x, l/2); // (1) in doc string
            }
            else if (m == l + 1)
            {
                val = x * (2*l + 1) * P_(l, l);  // (2)
            }
            else if (m <= l)
            {
                val = _consts_->get_leg_consts1_val(l, m) * x * P_(l-1, m) -
                            _consts_->get_leg_consts2_val(lInd, mInd) * P_(l-1, m);  // (3)
            }
            P_.set_val(lInd, mInd, val);
        }
    }
}


/*
 Calculate the spherical harmonics according to the equation:
 
    Y_(n,m)(theta, phi) = (-1)^m * sqrt((n-m)! / (n + m)!) * P_(n,m)(cos(theta)) * exp(i*m*phi)
 where P_(n, m) are the associated Legendre polynomials.
 
 */
void SHCalc::calc_sh()
{
    complex<double> iu (0, 1.0);  // complex unit
    int n, m;
    complex<double> val, mcomp;
    double shc;  // constant value
    for (n = 0; n < 2*nPoles_; n++)
    {
        for (m = 0; m < 2*nPoles_; m++)
        {
            shc = _consts_->get_sh_consts_val(n, m);
            mcomp = complex<double> (m, 0);
            val = pow(-1, (double) m) * shc * P_(n, m) * exp(iu * mcomp * phi_);
            Y_.set_val(n, m, val);
        }
    }
    
}

/*
 Return the results of the spherical harmonic calculation for an n, m.
 If m is negative, then we return the complex conjugate of the calculated 
 value for the positive value of m
 */
const complex<double> SHCalc::get_result(const int n, const int m) const
{
    if (m < 0)
    {
        return conj(Y_(n, m));  // complex conjugate
    }
    else
    {
        return Y_(n, m);
    }
}

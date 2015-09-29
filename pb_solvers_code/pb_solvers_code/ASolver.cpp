//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"


ASolver::ASolver(const vector<double>* a, const int N, const int p,
                const BesselCalc* _bcalc, const Constants* _consts)
:N_(N), p_(p), _besselCalc_(_bcalc), _consts_(_consts), _a_(a), gamma_(N, N)
,delta_(N, N), A_(N), E_(N)
{
}

/*
 Equation on Lotan 2006 page 544
 */
const double ASolver::calc_indi_gamma(int i, int n) const
{
    double g;  // output
    double kap = _consts_->get_kappa();
    double ai = _a_->at(i);  // radius
    double eps_p = _consts_->get_dielectric_prot();  // IS THIS WHAT epsilon_p is?? ##########################!!!!!!!!#####################
    double eps_s = _consts_->get_dielectric_water();
    vector<double> bk_all = _besselCalc_->calc_mbfK(n+1, kap*ai); // all bessel function k
    double bk2 = bk_all[n];   // bessel k at n+1
    double bk1 = bk_all[n-1];  // bessel k at n
    
    g = (2*n + 1) * exp(kap*ai);
    g = g / (((2*n + 1)*bk2) + (n*bk1*((eps_p / eps_s) - 1)));
    return g;
}

/*
 Equation on Lotan 2006 page 544
 */
const double ASolver::calc_indi_delta(int i, int n) const
{
    double d;  // output
    double kap = _consts_->get_kappa();
    double ai = _a_->at(i);  // radius
    double eps_p = _consts_->get_dielectric_prot();
    double eps_s = _consts_->get_dielectric_water();
    vector<double> bi_all = _besselCalc_->calc_mbfK(n+1, kap*ai); // all bessel function I
    double bi2 = bi_all[n];   // bessel i at n+1
    double bi1 = bi_all[n-1];  // bessel i at n
    
//    d = pow(ai, 2*n+1) / (2*n+1);
    d = (kap*kap*ai*ai) * (bi2 / (2*n + 3));
    d += (n * bi1 * (1 - (eps_p / eps_s)));
    d *= pow(ai, 2*n+1) / (2*n+1);
    return d;
}


/*
 Constructs the gamma matrix, which is a diagonal matrix of diagonal
 matrices that contain values calculated in calc_indi_gamma()
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
 matrices that contain values calculated in calc_indi_delta()
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

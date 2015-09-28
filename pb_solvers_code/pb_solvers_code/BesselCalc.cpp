//
//  BesselCalc.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "BesselCalc.h"


BesselConstants::BesselConstants(const int N)
:N_(N)
{
    recConsts_.reserve(2 * N_);
    int n;
    double val;
    for (n = 0; n < 2*N_; n++)
    {
        val = 1.0 / ((2 * n-1) * (2 * n-3));
        recConsts_.push_back(val);
    }
}


BesselCalc::BesselCalc(int N, BesselConstants* consts)
: N_(N), _consts_(consts)
{
    assert (_consts_->get_n() == N_);
}

const vector<double> BesselCalc::calc_mbfK(const int num_iter,
                                           const double z) const
{
    vector<double> K;
    K.reserve(num_iter);
    
    if (num_iter > 0) K.push_back(1.0);
    if (num_iter > 1) K.push_back(1 + z);
    
    double z_sq = z * z;
    int i;
    double val;
    for (i = 2; i < num_iter; i++)
    {
        val = K[i-1] + z_sq * K[i-2] * _consts_->get_const_val(i);
        K.push_back(val);
    }
    return K;
}

const vector<double> BesselCalc::calc_mbfI(const int num_iter,
                                           const double z) const
{
    vector<double> I;
    I.reserve(num_iter);
    for (int j = 0; j < num_iter; j++)
    {
        I.push_back(1);
    }
    
    if (z != 0)
    {
        double z2 = 0.5 * z * z;
        int k, m;
        double t;
        for (k = 0; k < num_iter; k++)
        {
            t = z2 / (2*k + 3);
            for (m = 0; m <= 20; m++)
            {
                I[k] += t;
                t *= z2 / ((m+1) * (2 * (k+m) + 3 ));  //EQ 1.15
                if (t < 1e-20) break;
            }
        }
    }
    return I;
}

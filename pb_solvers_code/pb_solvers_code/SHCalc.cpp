//
//  SHCalc.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "SHCalc.h"

SHCalc::SHCalc(int N)
:N_(N), consts1_(2*N, 2*N), consts2_(2*N, 2*N),
consts3_(2*N, 2*N), consts4_(2*N), IDX_(2*N_)
{
    vector<double> temp;
    temp.reserve(4 * N);
    temp.push_back(0);
    int i, n, m;
    for (i = 1; i < 4 * N_; i++)
    {
        temp.push_back(temp[i-1] * sqrt(i));
    }

    for (n = 0; n < 2 * N_; n++)
    {
        for (m = 0; m <= n; m++)
        {
            consts1_.set_val(n, m, (2*n-1) / (n-m));
            consts2_.set_val(n, m, (n+m-1) / (n-m));
            consts3_.set_val(n, m, temp[n-m] / temp[n+m]);
        }
    }
    
    consts4_[0] = 1.0;
    consts4_[1] = 1.0;
    IDX_[0] = 0;
    IDX_[1] = 1;
    for (i = 2; i < 2 * N_; i++)
    {
        consts4_[i] = consts4_[i-1] * (2*i - 1);
        IDX_[i] = IDX_[i] + i;
    }
    
}

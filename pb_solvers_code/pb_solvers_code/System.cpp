//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"


System::System(Constants consts, const vector<double>& a, const vector<double>& M,
               const vector<vector<double> > qs, const vector<EuPoint<double> > pos)
:a_(a), M_(M), qs_(qs), pos_(pos), consts_(consts)
{
}
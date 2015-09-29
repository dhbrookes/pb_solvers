//
//  System_hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef System_hpp
#define System_hpp

#include <stdio.h>
#include <vector>
#include "util.h"
#include "Constants.h"

using namespace std;
/*
 Class containing all of the relevant information for a particular system
 (molecule positions, charges, etc.)
 */

class System
{
protected:
    vector<double>               a_; // radii of every molecule
    vector<double>               M_; // number of charges in each molecule
    vector<vector<double> >      qs_; // magnitude of each charge in each molecule
    vector<EuPoint<double> >     pos_; // positions of every molecules
    Constants                    consts_;
    
public:
    System(Constants consts, const vector<double>& a, const vector<double>& M,
           const vector<vector<double> > qs, const vector<EuPoint<double> > pos);
    
    const double get_ai(int i) const                    { return a_[i];     }
    const double get_Mi(int i) const                    { return M_[i];     }
    const double get_qij(int i, int j)                  { return qs_[i][j]; }
    const EuPoint<double>& get_posi(int i) const        { return pos_[i];   }
    const Constants& get_consts() const                 { return consts_;   }
    const SphPoint<double> get_sph_posi(int i) const    { return pos_[i].convert_to_spherical(); }
    
};

#endif /* Setup_hpp */

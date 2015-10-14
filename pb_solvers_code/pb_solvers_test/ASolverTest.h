//
//  ASolverTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/8/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ASolverTest_h
#define ASolverTest_h

#include "ASolver.h"

class ASolverTest
{
public :
  void RunASolverTest()
  {
    Constants const_;
    vector< Molecule > mol_;
    
    mol_.clear( );
    EPt pos[2] = { EPt( 0.0, 0.0, 0.0 ), EPt( 0.0, 0.0, 5.0 ) };
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      double a = 1.0;
      vector<double> charges(1);
      vector<EPt> posCharges(1);
      
      charges[0] = 1.0;
      posCharges[0] = pos[molInd];
      
      Molecule molNew = Molecule( M, a, charges, posCharges );
      mol_.push_back( molNew );
    }
    
    const int vals           = nvals;
    BesselConstants bConsta  = BesselConstants( vals );
    BesselCalc bCalcu        = BesselCalc( vals, &bConsta);
    SHCalcConstants SHConsta = SHCalcConstants( vals );
    SHCalc SHCalcu           = SHCalc( vals, &SHConsta );
    System sys               = System( const_, mol_ );
    ASolver ASolvTest = ASolver( 2, vals, &bCalcu, &SHCalcu, sys);
    
  }
  
};



#endif /* ASolverTest_h */

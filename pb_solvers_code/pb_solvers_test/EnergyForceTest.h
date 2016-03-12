//
//  EnergyForceTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 1/8/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef EnergyForceTest_h
#define EnergyForceTest_h

#include "EnergyForce.h"

class EnForTest
{
  
  public :
  void RunEnForTest()
  {
    Constants const_;
    vector< Molecule > mol_;
    vector< Molecule > mol_sing_;
    
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.5 ), Pt( 11.0, 6.9, 24.3 ) };
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, 0.0 ) };
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0] = cg[molInd]; posCharges[0] = cgPos[molInd]; vdW[0] = 0.0;
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0; posCharges[0] = cgPosSi[molInd];
      Molecule molSing( "stat", 10.0, charges, posCharges, vdW);
      mol_sing_.push_back( molSing );
    }
    
    const int vals           = 10;
    BesselConstants bConsta( 2*vals );
    BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
    SHCalcConstants SHConsta( 2*vals );
    SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
    System sys( mol_ );
    ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                        sys.get_lambda(), vals);
    
    ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      make_shared<System> (sys),
                      make_shared<Constants> (const_), vals);
    ASolvTest.solve_A( 1E-12 );
    ASolvTest.solve_gradA(1E-12);

    ForceCalc FoTest( ASolvTest );
    TorqueCalc TorTest( ASolvTest );
    
  }

};


#endif /* EnergyForceTest_h */

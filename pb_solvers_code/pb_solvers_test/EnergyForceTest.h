//
//  EnergyForceTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 1/8/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef EnergyForceTest_h
#define EnergyForceTest_h

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
      vector<double> charges(1);
      vector<Pt> posCharges(1);
      
      charges[0]    = cg[molInd];
      posCharges[0] = cgPos[molInd];
      
      Molecule molNew( M, rd[molInd], charges, posCharges, pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0;
      posCharges[0] = cgPosSi[molInd];
      
      Molecule molSing( M, 10.0, charges, posCharges);
      mol_sing_.push_back( molSing );
    }
    
    const int vals           = 10;
    BesselConstants bConsta( 2*vals );
    BesselCalc bCalcu( 2*vals, bConsta);
    SHCalcConstants SHConsta( 2*vals );
    SHCalc SHCalcu( 2*vals, SHConsta );
    System sys( const_, mol_ );
    ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                        sys.get_lambda(), nvals);
    
    ASolver ASolvTest( 2, vals, bCalcu, SHCalcu, sys);
    ASolvTest.solve_A( 1E-12 );
    ASolvTest.solve_gradA(1E-12);
    
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, 2, nvals);
    
//    for (int n=0; n<2; n++)
//      cout << " This is my force " << n << " : " << FoTest.get_fi(n)[0] << ", " << FoTest.get_fi(n)[1] << ", " << FoTest.get_fi(n)[2] << endl;
    
  }
  
};


#endif /* EnergyForceTest_h */

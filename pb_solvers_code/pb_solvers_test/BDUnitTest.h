//
//  BDUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 1/25/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef BDUnitTest_h
#define BDUnitTest_h

#include "BD.h"
#include "EnergyForce.h"

class BDUTest : public ::testing::Test
{
public :
  
protected :
  int vals_;
  Constants const_;
  vector< Molecule > mol3_; vector< Molecule > mol_;
  vector< Molecule > mol_sing_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    double cg[2]  = { 5.0, -0.4}; double rd[2] = { 5.6, 10.4};
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, 0.0 ) };
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(1); vector<Pt> posCharges(1);
      charges[0]    = cg[molInd]; posCharges[0] = cgPos[molInd];
      
      Molecule molNew( M, rd[molInd], charges, posCharges, pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0; posCharges[0] = cgPosSi[molInd];
      Molecule molSing( M, 10.0, charges, posCharges);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  virtual void TearDown() {}
} ; // end BDUTest

TEST_F(BDUTest, Move)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, -5.0), Pt(10.0, 7.8, 25.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 1; vector<double> charges(1); vector<Pt> posCharges(1);
    charges[0] = 2.0; posCharges[0] = pos[molInd];
    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
    mol_.push_back( molNew );
  }
  const int vals           = nvals;
  int nmol                 = 3;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
                   ASolvTest.calc_L(), ASolvTest.calc_gradL(),
                   const_, nmol, nvals);
  
  vector<double> dTr(nmol); vector<double> dRot(nmol);
  dTr[0]  = 0.93; dTr[1]  = 0.93; dTr[2]  = 0.93;
  dRot[0] = 0.93; dRot[1] = 0.93; dRot[2] = 0.93;
  
  BD BDTest( sys, 2.0, dTr, dRot);
}

#endif /* BDUnitTest_h */

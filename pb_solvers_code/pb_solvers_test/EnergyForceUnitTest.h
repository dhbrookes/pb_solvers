//
//  EnergyForceUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 1/7/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef EnergyForceUnitTest_h
#define EnergyForceUnitTest_h

#include "EnergyForce.h"

class EnergyForceUTest : public ::testing::Test
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
    double cg[2] = { 5.0, -0.4}; double rd[2] = { 5.6, 10.4};
    
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
  } // end SetUp
  
  virtual void TearDown() {}
  
  double Mol1F[3] = {4.13649103e-17, -1.72304216e-05, -6.62708522e-05};
  double Mol2F[3] = {0.000124230113, 8.60763311e-06, 3.31062812e-05};
  double Mol3F[3] = {-0.000124230111, 8.60763311e-06, 3.31062812e-05};

  
} ; // end EnergyForceUTest

TEST_F(EnergyForceUTest, checkEnergy)
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
  ASolvTest.solve_A(1E-20);
  
  EnergyCalc EnTest( ASolvTest.get_A(), ASolvTest.calc_L(),
                     const_, nmol, nvals);

  EXPECT_NEAR( EnTest.get_omega_i_int(0), 0.00117849131, preclim);
  EXPECT_NEAR( EnTest.get_omega_i_int(1), 0.00199508625, preclim);
  EXPECT_NEAR( EnTest.get_omega_i_int(2), 0.00199508625, preclim);
}

TEST_F(EnergyForceUTest, checkEnergySing)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20);
  
  EnergyCalc EnTest( ASolvTest.get_A(), ASolvTest.calc_L(),
                    const_, nmol, nvals);
  
  for (int n=0; n<2; n++)
    EXPECT_NEAR( EnTest.get_omega_i_int(0), 0.000573165, preclim);
}


TEST_F(EnergyForceUTest, checkForce)
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
  
  for (int n=0; n<3; n++)
  {
    EXPECT_NEAR( FoTest.get_fi(0)[n]/Mol1F[n], 1.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(1)[n]/Mol2F[n], 1.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(2)[n]/Mol3F[n], 1.0, preclim);
  }
}

TEST_F(EnergyForceUTest, checkForceSing)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
                   ASolvTest.calc_L(), ASolvTest.calc_gradL(),
                   const_, nmol, nvals);
  
  for (int n=0; n<2; n++)
  {
    EXPECT_NEAR( FoTest.get_fi(n)[0], 0.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(n)[1], 0.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(n)[2], -3.60903749e-05*pow(-1.0,n), preclim);
  }
}

  
#endif /* EnergyForceUnitTest_h */

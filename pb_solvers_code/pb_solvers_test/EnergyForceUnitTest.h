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
      int M = 1; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      
      Molecule molNew( "stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]=2.0; vdW[0] = 0.0; posCharges[0] = cgPosSi[molInd];
      Molecule molSing( "stat", 10.0, charges, posCharges, vdW);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  double MolTripSing[3] = {0.165835788, 0.11870012, 0.11870012};
  
  double Mol1F[3] = {-1.11207811e-07, -0.000154561614, -0.000595623838};
  double Mol2F[3] = {0.00112021774, 8.12434314e-05, 0.00029752783};
  double Mol3F[3] = {-0.00112010937, 7.32743533e-05, 0.000297943745};
  
  double Tor1[3] = {-0.000203298157, 0.000195890873, -5.18379606e-05};
  double Tor2[3] = {9.82164529e-05, -9.71696398e-05, -0.000371317558};
  double Tor3[3] = {9.52340348e-05, -0.00010371257, 0.000369189207};
 
} ; // end EnergyForceUTest

TEST_F(EnergyForceUTest, checkEnergy)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, -5.0), Pt(10.0, 7.8, 25.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 1; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; vdW[0] = 0.0; posCharges[0] = pos[molInd];
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  const int vals           = nvals;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  shared_ptr<System> sys = make_shared<System>( mol_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys->get_lambda(), vals);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    sys,
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-20);
  
  EnergyCalc EnTest( make_shared<ASolver> (ASolvTest));

  EXPECT_NEAR( EnTest.get_omega_i_int(0), 0.00117849131, preclim);
  EXPECT_NEAR( EnTest.get_omega_i_int(1), 0.00199508625, preclim);
  EXPECT_NEAR( EnTest.get_omega_i_int(2), 0.00199508625, preclim);
}

TEST_F(EnergyForceUTest, checkEnergySing)
{
  const int vals           = nvals;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    make_shared<System> (sys),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-20);
  
  EnergyCalc EnTest( make_shared<ASolver> (ASolvTest));
  
  for (int n=0; n<mol_sing_.size(); n++)
    EXPECT_NEAR( EnTest.get_omega_i_int(0), 0.000573165, preclim);
}

TEST_F(EnergyForceUTest, checkEnergySingMulti)
{
  mol_sing_.clear( );
  Pt pos[3] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 0.0, 0.0, -5.0 ),Pt( 0.0, 0.0, 5.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_sing_.push_back( molNew );
  }
  const int vals           = nvals;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    make_shared<System> (sys),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-20);
  
  EnergyCalc EnTest( make_shared<ASolver> (ASolvTest) );
  
  for (int n=0; n<mol_.size(); n++)
    EXPECT_NEAR( EnTest.get_omega_i_int(n)/MolTripSing[n], 1, preclim);
}

TEST_F(EnergyForceUTest, checkForce)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, -5.0), Pt(10.0, 7.8, 25.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  const int vals           = nvals;
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
  ASolvTest.solve_A(1E-40); ASolvTest.solve_gradA(1E-40);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest) );
  
  for (int n=0; n<mol_.size(); n++)
  {
    EXPECT_NEAR( FoTest.get_fi(0)[n]/Mol1F[n], 1.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(1)[n]/Mol2F[n], 1.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(2)[n]/Mol3F[n], 1.0, preclim);
  }
}

TEST_F(EnergyForceUTest, checkForceSing)
{
  const int vals           = nvals;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    make_shared<System> (sys),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-40); ASolvTest.solve_gradA(1E-40);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest) );

  for (int n=0; n<2; n++)
  {
    EXPECT_NEAR( FoTest.get_fi(n)[0], 0.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(n)[1], 0.0, preclim);
    EXPECT_NEAR( FoTest.get_fi(n)[2]/(-3.61229952e-05*pow(-1.0,n)), 1, preclim);
  }
}

TEST_F(EnergyForceUTest, checkForce3Cg)
{
  mol_.clear( );
  Pt pos[2] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 0.0, -5.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  const int vals           = 5;
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
  ASolvTest.solve_A(1E-30); ASolvTest.solve_gradA(1E-30);

  EnergyCalc EnTest( make_shared<ASolver> (ASolvTest));
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  
  EXPECT_NEAR( EnTest.get_omega_i_int(0)/0.0845178625, 1, preclim);
  EXPECT_NEAR( FoTest.get_fi(0)[0], 0, 1e-12);
  EXPECT_NEAR( FoTest.get_fi(0)[1], 0, 1e-12);
  EXPECT_NEAR( FoTest.get_fi(0)[2]/0.0239714121, 1.0, preclim);
  EXPECT_NEAR( EnTest.get_omega_i_int(1)/0.0845178625, 1, preclim);
  EXPECT_NEAR( FoTest.get_fi(1)[0], 0, 1e-12);
  EXPECT_NEAR( FoTest.get_fi(1)[1], 0, 1e-12);
  EXPECT_NEAR( FoTest.get_fi(1)[2]/-0.0239714121, 1.0, preclim);

}

TEST_F(EnergyForceUTest, checkTorque)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, -5.0), Pt(10.0, 7.8, 25.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  const int vals           = 5;
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
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
  
  for (int n=0; n<3; n++)
  {
    EXPECT_NEAR( TorTest.get_taui(0)[n]/Tor1[n], 1, preclim);
    EXPECT_NEAR( TorTest.get_taui(1)[n]/Tor2[n], 1, preclim);
    EXPECT_NEAR( TorTest.get_taui(2)[n]/Tor3[n], 1, preclim);
  }
}

TEST_F(EnergyForceUTest, checkTorqueSing)
{
  const int vals           = nvals;
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
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
  
  for (int n=0; n<mol_.size(); n++)
    for (int dim=0; dim<3; dim++)
      EXPECT_NEAR( TorTest.get_taui(n)[dim], 0.0, 1e-14);
}

TEST_F(EnergyForceUTest, checkTorqueSing3)
{
  mol_sing_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0),Pt(0.0, 0.0, -5.0),Pt(0.0, 0.0, 5.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_sing_.push_back( molNew );
  }

  const int vals           = 5;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    make_shared<System> (sys),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));

  EXPECT_NEAR( TorTest.get_taui(0)[0], 0.0, 1e-14);
  EXPECT_NEAR( TorTest.get_taui(0)[1], 0.0, 1e-14);
  EXPECT_NEAR( TorTest.get_taui(0)[2], 0.0, 1e-14);
  
  EXPECT_NEAR( TorTest.get_taui(1)[0], -0.00799959051, 1e-9);
  EXPECT_NEAR( TorTest.get_taui(1)[1],  0.00799959051, 1e-9);
  EXPECT_NEAR( TorTest.get_taui(1)[2], 0.0, 1e-14);
  
  EXPECT_NEAR( TorTest.get_taui(2)[0],  0.00799959051, 1e-9);
  EXPECT_NEAR( TorTest.get_taui(2)[1], -0.00799959051, 1e-9);
  EXPECT_NEAR( TorTest.get_taui(2)[2], 0.0, 1e-14);

}

  
#endif /* EnergyForceUnitTest_h */

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
      vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0] = cg[molInd]; posCharges[0] = cgPos[molInd]; vdW[0] = 0.0;
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0; posCharges[0] = cgPosSi[molInd];
      Molecule molSing( "stat", 10.0, charges, posCharges, vdW);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  double BD1Force[11] = {5.01987155,5.03957651,5.05911838,5.07850055,5.0977263,
    5.11679881,5.13572112,5.15449623,5.17312699,5.19161619,5.20996653};
  
  double BD1NegFo[11] = {4.9813721,4.96261031,4.94371229,4.92467563,4.90549784,
    4.88617636,4.86670858,4.84709179,4.8273232,4.80739997,4.78731914};
  
  virtual void TearDown() {}
} ; // end BDUTest

// Test random vector? test check for collision, trans update, rot update

TEST_F(BDUTest, ForcePos)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd]);
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
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dTr[0]  = 0.0;  dTr[1]  = 0.01;
  dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
            dTr, dRot, false);
  
  for (int step=0; step < 10; step ++)
  {
    ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      BDTest.get_system(),
                      make_shared<Constants> (const_), vals);
    ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
    
    ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
    TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, nmol, vals);
//    TorqueCalc TorTest( SHCalcu, bCalcu, ASolvTest.calc_gradL(),
//                       ASolvTest.get_gamma(), const_, sys, vals);
    
    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
    
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).y()
                /BD1Force[step], 1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
  }
}

TEST_F(BDUTest, ForcePosZ)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 0.0, 5.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
//    int M = 1; vector<double> charges(1); vector<Pt> posCharges(1);
//    charges[0] = 2.0; posCharges[0] = pos[molInd];
//    Molecule molNew( M, 1.0, charges, posCharges, pos[molInd]);
//    mol_.push_back( molNew );
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd]);
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
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dTr[0]  = 0.0;  dTr[1]  = 0.01;
  dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
            dTr, dRot, false);
  
  for (int step=0; step < 10; step ++)
  {
    ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      BDTest.get_system(),
                      make_shared<Constants> (const_), vals);
    ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
    
    ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
    TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, nmol, vals);
//    TorqueCalc TorTest( SHCalcu, bCalcu, ASolvTest.calc_gradL(),
//                       ASolvTest.get_gamma(), const_, sys, vals);
    
    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
    
    
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).z()
                /BD1Force[step], 1.0, preclim);
  }
}

TEST_F(BDUTest, ForceOpp)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0*pow(-1, molInd); posCharges[0]=pos[molInd]; vdW[0]=0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd]);
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
  
  
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dTr[0]  = 0.0;  dTr[1]  = 0.01;
  dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
            dTr, dRot, false);
  
  for (int step=0; step < 10; step ++)
  {
    ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      BDTest.get_system(), 
                      make_shared<Constants> (const_), vals);
    ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
    
    ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
    TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
    
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, nmol, vals);
//    TorqueCalc TorTest( SHCalcu, bCalcu, ASolvTest.calc_gradL(),
//                       ASolvTest.get_gamma(), const_, sys, vals);
    
    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
    
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(1).y()
                /BD1NegFo[step], 1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
  }
}

TEST_F(BDUTest, TorquePos)
{
  mol_.clear( );
  Pt pos[2] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 3;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd];
    charges[1]=2.0; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
    charges[2]=2.0; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);
    Molecule molNew( "stat", 1.4714, charges, posCharges, vdW);
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
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dTr[0]  = 0.01;  dTr[1]  = 0.01;
  dRot[0] = 0.01; dRot[1] = 0.01;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
            dTr, dRot, false);
  
  for (int step=0; step < 10; step ++)
  {
    ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      make_shared<System> (sys),
                      make_shared<Constants> (const_), vals);
    ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
    
    ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
    TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, nmol, vals);
//    
//    TorqueCalc TorTest( SHCalcu, bCalcu, ASolvTest.calc_gradL(),
//                       ASolvTest.get_gamma(), const_, sys, vals);
  
    double f0 = (abs(TorTest.get_taui(0)[0])<1e-15) ? 0:TorTest.get_taui(0)[0];
    double f1 = (abs(TorTest.get_taui(0)[1])<1e-15) ? 0:TorTest.get_taui(0)[1];
    double f2 = (abs(TorTest.get_taui(0)[2])<1e-15) ? 0:TorTest.get_taui(0)[2];
    cout << " This is tau " << f0 << " " << f1 << " " << f2 << " " << endl;
    f0 = (abs(TorTest.get_taui(1)[0])<1e-15) ? 0:TorTest.get_taui(1)[0];
    f1 = (abs(TorTest.get_taui(1)[1])<1e-15) ? 0:TorTest.get_taui(1)[1];
    f2 = (abs(TorTest.get_taui(1)[2])<1e-15) ? 0:TorTest.get_taui(1)[2];
    cout << " This is tau " << f0 << " " << f1 << " " << f2 << " " << endl;
    
    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
  
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).x(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).y(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).z(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(1).x(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(1).y()
//                /BD1NegFo[step], 1.0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).z(), 0, preclim);
  }
}


//TEST_F(BDUTest, TorqueOpp)
//{
//  mol_.clear( );
//  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0), Pt(-10.0, 7.8, 25.0)};
//  for (int molInd = 0; molInd < 2; molInd ++ )
//  {
//    int M = 3; vector<double> charges(M); vector<Pt> posCharges(M);
//    charges[0] = 2.0*pow(-1, molInd); posCharges[0] = pos[molInd];
//    charges[1] = 2.0*pow(-1, molInd); posCharges[1] = pos[molInd]+Pt(1,0,0);
//    charges[2] = 2.0*pow(-1, molInd); posCharges[2] = pos[molInd]+Pt(0,1,0);
//    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
//    mol_.push_back( molNew );
//  }
//  const int vals           = 5;
//  int nmol                 = 2;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), vals);
//  
//  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
//  dTr[0]  = 0.0;  dTr[1]  = 0.01;
//  dRot[0] = 0.0; dRot[1] = 0.01;
//  BD BDTest( sys, dTr, dRot, false);
//  
//  for (int step=0; step < 10; step ++)
//  {
//    ASolver ASolvTest( nmol, vals, bCalcu, SHCalcu, BDTest.get_system());
//    ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
//    
//    ForceCalc FoTest( ASolvTest.get_A(), ASolvTest.get_gradA(),
//                     ASolvTest.calc_L(), ASolvTest.calc_gradL(),
//                     const_, nmol, vals);
//    
//    TorqueCalc TorTest( SHCalcu, bCalcu, ASolvTest.calc_gradL(),
//                       ASolvTest.get_gamma(), const_, sys, vals);
//    
//    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
//    
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).x(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).y(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).z(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(1).x(), 0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(1).y()
//                /BD1NegFo[step], 1.0, preclim);
//    EXPECT_NEAR(BDTest.get_system().get_centeri(0).z(), 0, preclim);
//  }
//}

#endif /* BDUnitTest_h */

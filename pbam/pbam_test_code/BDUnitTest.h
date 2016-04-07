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
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd],
                      molInd, 0);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0; posCharges[0] = cgPosSi[molInd];
      Molecule molSing( "stat", 10.0, charges, posCharges, vdW, molInd, 0);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  double BD1Force[11] = {5.01987155,5.03957651,5.05911838,5.07850055,5.0977263,
    5.11679881,5.13572112,5.15449623,5.17312699,5.19161619,5.20996653};
  double BD1NegFo[11] = {4.9813721,4.96261031,4.94371229,4.92467563,4.90549784,
    4.88617636,4.86670858,4.84709179,4.8273232,4.80739997,4.78731914};
  
  double BDTor0x[11] = {0.511283938,0.521770587,0.531550359,0.540697963,
    0.549276092,0.557338043,0.564929619,0.572090556,0.578855605,0.585255382};
  double BDTor0y[11] = {0.488455458,0.477237315,0.466319865,0.455681592,
    0.445304137,0.435171582,0.425269944,0.415586809,0.406111054,0.396832632};
  double BDTor1x[11] = {-0.511283938,-0.521770587,-0.531550359,-0.540697963,
    -0.54927609,-0.55733804,-0.56492962,-0.57209056,-0.57885561,-0.58525538};
  double BDTor1y[11] = {-0.488455458,-0.477237315,-0.466319865,-0.455681592,
    -0.445304137,-0.43517158,-0.42526994,-0.41558681,-0.40611105,-0.39683263};
  
  double BD3Tor01x[11] = {0.999737677, 0.99899179, 0.997819937, 0.996274456,
    0.994402575,0.992246596,0.989844074,0.987227999,0.984426951,0.981465243};
  double BD3Tor01y[11] = {0.0228681166,0.0448222486,0.0658891503,0.0860987846,
    0.105483774,0.124078987,0.141921255,0.159049193,0.175503109,0.191324943};
  double BD3Tor01z[11] = {-0.00127504644,-0.00252396704,-0.0037407794,
    -0.0049201464,-0.00605734744,-0.00714826235,-0.00818937227,
    -0.0091777811,-0.0101112586,-0.0109883033};
  double BD3Tor02x[11] = {-0.0228699223,-0.0448290536,-0.0659035515,
    -0.0861228226,-0.105518973,-0.124126394,-0.141981481,-0.159122453,
    -0.175589265,-0.191423544};
  double BD3Tor02y[11] = {0.999737441,0.998991004,0.997818506,0.996272466,
    0.994400266,0.992244333,0.989842323,0.987227306,0.984427922,0.981468526};
  double BD3Tor02z[11] = {-0.00142003996,-0.0027073701,-0.00386660605,
    -0.00490224967,-0.00581863072,-0.00661987474,-0.00730990146,-0.0078924592,
    -0.00837120082,-0.00874980532};
  double BD3Tor11x[11] = {0.998988091,0.995741444,0.989923972,0.981176582,
    0.969123764,0.953383976,0.933584726,0.909383147,0.880492404,0.846713419};
  double BD3Tor11y[11] = {-0.0449536544,-0.0921480821,-0.141540056,-0.193036939,
    -0.246485329,-0.301659001,-0.358247005,-0.41584335,-0.473940303,-0.531928025};
  double BD3Tor11z[11] = {0.0014009485,0.00277625458,0.00411601907,
    0.00540882609,0.00664169132,0.0078000842,0.00886805976,0.00982854209,
    0.010663804,0.0113561803};
  double BD3Tor12x[11] = {0.0449514378,0.0921390033,0.141519152,0.192998941,
    0.246424682,0.301569898,0.358123443,0.415679198,0.473729397,0.531664279};
  double BD3Tor12y[11] = {0.998987931,0.99574089,0.989922933,0.981175119,
    0.969122104,0.953382514,0.933584024,0.909383926,0.880495521,0.846719835};
  double BD3Tor12z[11] = {0.00157546932,0.00323784502,0.00499163379,
    0.00684065397,0.00878770676,0.0108341784,0.0129795761,0.0152210136,
    0.017552679,0.0199653419};
  double BD3Tor21x[11] = {0.999999953,0.999999816,0.999999591,0.999999281,
    0.99999889,0.999998422,0.999997878,0.999997264,0.999996582,0.999995834};
  double BD3Tor21y[11] = {0.000248104461,0.000493698151,0.000736689366,
    0.000976983779,0.00121448613,0.00144910249,0.00168074309,0.00190932586,
    0.00213478061,0.00235705375};
  double BD3Tor21z[11] = {-0.000177584848,-0.000352645416,-0.000525193586,
    -0.000695250272,-0.000862847668,-0.00102803178,-0.00119086513,
    -0.00135142956,-0.00150982881,-0.00166619074};
  double BD3Tor22x[11] = {-0.000248141178,-0.00049384295,-0.000737010545,
    -0.000977546636,-0.00121535305,-0.00145033306,-0.00168239421,-0.00191145198,
    -0.00213743388,-0.00236028431};
  double BD3Tor22y[11] = {0.999999948,0.999999794,0.999999541,0.999999194,
    0.999998757,0.999998232,0.999997623,0.999996935,0.999996171,0.999995334};
  double BD3Tor22z[11] = {-0.00020676628,-0.0004106376,-0.000611612002,
    -0.000809695912,-0.00100490659,-0.0011972749,-0.00138684841,-0.0015736946,
    -0.00175790403,-0.00193959316};
  
  double bdRun10[2] = {-3.77575521, 9.77575521};
  
  double bdRun30x[2] = {-0.634593964,0.634593964};
  double bdRun30y[2] = {-20.8011211,26.8011211};

  double bdRun30Rot[2][3][3] = {{{0,0,0},{-0.800362355,0.565300727,0.199637645},{0.199637645, 0.565300727, -0.800362355}}, {{0,0,0},
      {-0.800362773, -0.565300282, 0.199637227},{0.199637227, -0.565300282,
        -0.800362773}}};
  
  virtual void TearDown() {}
} ; // end BDUTest

// Test random vector? test check for collision, trans update, rot update
TEST_F(BDUTest, dtTest)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.0;  dTr[1]  = 0.01; dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
                dTr, dRot, false);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    BDTest.get_system(),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-4); ASolvTest.solve_gradA(1E-4);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
  FoTest.calc_force(); TorTest.calc_tau();
  
  BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
  
  EXPECT_NEAR(BDTest.get_dt()/2.0, 1.0, preclim);
}

TEST_F(BDUTest, dtLargeTest)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.0;  dTr[1]  = 0.01; dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
                dTr, dRot, false);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    BDTest.get_system(),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-4); ASolvTest.solve_gradA(1E-4);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
  FoTest.calc_force(); TorTest.calc_tau();
  
  BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
  
  EXPECT_NEAR(BDTest.get_dt()/2.535522526, 1.0, preclim);
}

TEST_F(BDUTest, distPBCTest)
{
  mol_.clear( );
  Pt pos[2] = {Pt(0.0, 0.0, 0.0), Pt(-25.0, 55.8, 21.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }
  double cutoff  = 10.00;
  double boxl    = 20.00;
  const int vals = 5;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_, cutoff, boxl );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dTr[0]  = 0.0;  dTr[1]  = 0.01; dRot[0] = 0.0; dRot[1] = 0.0;
  BDStep BDTest( make_shared<System> (sys), make_shared<Constants> (const_),
                dTr, dRot, false);
  
  ASolver ASolvTest(make_shared<BesselCalc> (bCalcu),
                    make_shared<SHCalc> (SHCalcu),
                    BDTest.get_system(),
                    make_shared<Constants> (const_), vals);
  ASolvTest.solve_A(1E-4); ASolvTest.solve_gradA(1E-4);
  
  ForceCalc FoTest( make_shared<ASolver> (ASolvTest));
  TorqueCalc TorTest( make_shared<ASolver> (ASolvTest));
  FoTest.calc_force(); TorTest.calc_tau();
  
  BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
  
  EXPECT_NEAR(BDTest.get_dt()/2.0, 1.0, preclim);
}


TEST_F(BDUTest, ForcePos)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.0;  dTr[1]  = 0.01; dRot[0] = 0.0; dRot[1] = 0.0;
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
    FoTest.calc_force(); TorTest.calc_tau();
    
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
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.0;  dTr[1]  = 0.01; dRot[0] = 0.0; dRot[1] = 0.0;
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
    FoTest.calc_force(); TorTest.calc_tau();
    
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
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.0;  dTr[1]  = 0.01;   dRot[0] = 0.0; dRot[1] = 0.0;
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

    FoTest.calc_force();  TorTest.calc_tau();
    
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
  Pt chgLoc[2] = {Pt(0.5, 0.5, 0.0), Pt(-0.5, -0.5, 0.0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0]=2.0; vdW[0]=0.0; posCharges[0]=pos[molInd]+chgLoc[molInd];
    Molecule molNew( "rot", 1.4714, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
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
  dTr[0]  = 0.01;  dTr[1]  = 0.01; dRot[0] = 0.01; dRot[1] = 0.01;
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
    FoTest.calc_force(); TorTest.calc_tau();

    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());
  
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).x()/BDTor0x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).y()/BDTor0y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).z(), 0, preclim);

    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).x()/BDTor1x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).y()/BDTor1y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).z(), 0, preclim);
  }
}


TEST_F(BDUTest, TorqueOpp)
{
  mol_.clear( );
  Pt pos[3] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0), Pt(-10.0, 7.8, 25.0)};
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0*pow(-1, molInd); posCharges[0] = pos[molInd];
    charges[1] = 2.0*pow(-1, molInd); posCharges[1] = pos[molInd]+Pt(1,0,0);
    charges[2] = 2.0*pow(-1, molInd); posCharges[2] = pos[molInd]+Pt(0,1,0);
    vdW[0]=0.0; vdW[1]=0.0; vdW[2]=0.0;
    Molecule molNew( "rot", 2.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }
  const int vals = 5;
  BesselConstants bConsta( 2*vals );
  BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
  SHCalcConstants SHConsta( 2*vals );
  SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
  System sys( mol_ );
  ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                      sys.get_lambda(), vals);
  
  vector<double> dTr(mol_.size()); vector<double> dRot(mol_.size());
  dRot[0] = 0.01; dRot[1] = 0.01; dRot[2] = 0.01;
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
    FoTest.calc_force(); TorTest.calc_tau();
    
    BDTest.bd_update(FoTest.get_F(), TorTest.get_Tau());

    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_centeri(0).z(), 0, preclim);
    
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 1).x()/BD3Tor01x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 1).y()/BD3Tor01y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 1).z()/BD3Tor01z[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 2).x()/BD3Tor02x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 2).y()/BD3Tor02y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(0, 2).z()/BD3Tor02z[step],
                1.0, preclim);
    
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 1).x()/BD3Tor11x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 1).y()/BD3Tor11y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 1).z()/BD3Tor11z[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 2).x()/BD3Tor12x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 2).y()/BD3Tor12y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(1, 2).z()/BD3Tor12z[step],
                1.0, preclim);
    
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 0).x(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 0).y(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 0).z(), 0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 1).x()/BD3Tor21x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 1).y()/BD3Tor21y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 1).z()/BD3Tor21z[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 2).x()/BD3Tor22x[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 2).y()/BD3Tor22y[step],
                1.0, preclim);
    EXPECT_NEAR(BDTest.get_system()->get_posij(2, 2).z()/BD3Tor22z[step],
                1.0, preclim);
  }
}

TEST_F(BDUTest, BDrunTimeTermY)
{
  vector<Molecule> mol;
  const int ml = 2;
  Pt pos[ml] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0)};
  for (int mi = 0; mi < ml; mi ++ )
  {
    int M = 1; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 7.0; posCharges[0] = pos[mi]; vdW[0]=0.0;
    Molecule molNew( "trans", 2.0, charges, posCharges, vdW, pos[mi],
                    mi, 0, 0, 0.1);
    mol.push_back( molNew );
  }
  const int vals = 5;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  auto ASolvTest = make_shared<ASolver>( bCalcu, SHCalcu, sys,
                                        make_shared<Constants> (const_), vals);

  shared_ptr<TimeTerminate> term = make_shared<TimeTerminate>(10);
  BDRun BDTest( ASolvTest, term, "", 0, false, true, 1e7, 1e-20);
  BDTest.run();
  
  EXPECT_NEAR(sys->get_time()/10, 1, preclim);
  for (int mi = 0; mi < ml; mi ++ )
  {
    EXPECT_NEAR(sys->get_centeri(mi).x(), 0, preclim);
    EXPECT_NEAR(sys->get_centeri(mi).y()/bdRun10[mi], 1, preclim);
    EXPECT_NEAR(sys->get_centeri(mi).z(), 0, preclim);
  }
}

TEST_F(BDUTest, BDrunTimeTermXY)
{
  vector<Molecule> mol;
  const int ml = 2;
  Pt pos[ml] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0)};
  for (int mi = 0; mi < ml; mi ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 7.0; posCharges[0] = pos[mi]; vdW[0]=0.0;
    charges[1] = 7.0; posCharges[1] = pos[mi]+Pt(1,0,0); vdW[1]=0.0;
    charges[2] = 7.0; posCharges[2] = pos[mi]+Pt(0,1,0); vdW[2]=0.0;
    Molecule molNew( "trans", 2.0, charges, posCharges, vdW, pos[mi],
                    mi, 0, 0, 0.1);
    mol.push_back( molNew );
  }
  const int vals = 5;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  auto ASolvTest = make_shared<ASolver>( bCalcu, SHCalcu, sys,
                                  make_shared<Constants> (const_), vals);
  
  shared_ptr<TimeTerminate> term = make_shared<TimeTerminate>(30);
  BDRun BDTest( ASolvTest, term, "", 0, false, true, 1e7, 1e-20);
  BDTest.run();
  
  EXPECT_NEAR(sys->get_time()/31.361344, 1, preclim);
  
  for (int mi = 0; mi < ml; mi ++ )
  {
    EXPECT_NEAR(sys->get_centeri(mi).x()/bdRun30x[mi], 1, preclim);
    EXPECT_NEAR(sys->get_centeri(mi).y()/bdRun30y[mi], 1, preclim);
    EXPECT_NEAR(sys->get_centeri(mi).z(), 0, preclim);
  }
}

TEST_F(BDUTest, BDrunTimeTermRot)
{
  vector<Molecule> mol;
  const int ml = 2;
  Pt pos[ml] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0)};
  for (int mi = 0; mi < ml; mi ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 7.0; posCharges[0] = pos[mi]; vdW[0]=0.0;
    charges[1] = 7.0; posCharges[1] = pos[mi]+Pt(0,0,1); vdW[1]=0.0;
    charges[2] = 7.0; posCharges[2] = pos[mi]+Pt(1,0,0); vdW[2]=0.0;
    Molecule molNew( "rot", 2.0, charges, posCharges, vdW, pos[mi], mi, 0,
                    0.1, 0.0);
    mol.push_back( molNew );
  }
  const int vals = 5;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  auto ASolvTest = make_shared<ASolver>( bCalcu, SHCalcu, sys,
                                        make_shared<Constants> (const_), vals);
  
  shared_ptr<TimeTerminate> term = make_shared<TimeTerminate>(30);
  BDRun BDTest( ASolvTest, term, "", 0, false, true, 1e7, 1e-30);
  BDTest.run();
  
  for (int mi = 0; mi < ml; mi ++ )
  {
    EXPECT_NEAR(sys->get_centeri(mi).x(), 0, preclim);
    if ( pos[mi].y() != 0)
      EXPECT_NEAR(sys->get_centeri(mi).y()/pos[mi].y(), 1, preclim);
    EXPECT_NEAR(sys->get_centeri(mi).z(), 0, preclim);
    
    for (int j=0; j < sys->get_Mi(mi); j++)
    {
      if (bdRun30Rot[mi][j][0] != 0)
        EXPECT_NEAR(sys->get_posij(mi,j).x()/bdRun30Rot[mi][j][0],1,preclim);
      if (bdRun30Rot[mi][j][1] != 0)
        EXPECT_NEAR(sys->get_posij(mi,j).y()/bdRun30Rot[mi][j][1],1,preclim);
      if (bdRun30Rot[mi][j][2] != 0)
        EXPECT_NEAR(sys->get_posij(mi,j).z()/bdRun30Rot[mi][j][2],1,preclim);
    }
  }
}

#endif /* BDUnitTest_h */

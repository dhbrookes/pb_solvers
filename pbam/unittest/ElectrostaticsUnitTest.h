//
//  ElectrostaticsUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ElectrostaticsUnitTest_h
#define ElectrostaticsUnitTest_h

#include <stdio.h>
#include "Electrostatics.h"

class ElecUTest : public ::testing::Test
{
public :

protected :
  
  int vals_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -8.0 ), Pt( 0, 0, 0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -7.0 ), Pt( 0, 0, 1 ) };
    double cg[2] = { -5.0, -5.0};
    double rd[2] = {  3.6,  3.6};
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 3; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      charges[1]=cg[molInd]; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
      charges[2]=cg[molInd]; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);

      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd],
                      molInd, 0);
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  double min[3] = {-7.2, -7.2, -15.2};
  double bin[3] = {1.30909091, 1.30909091, 2.03636364};
  double bn2[3] = {0.72, 0.72, 1.12};
  
  double potE[300] = {-0.0152695428,-0.0173031116,-0.0193525353,-0.0212205,-0.0226729181,-0.0235165772,-0.0236615797,-0.0231151987,-0.0219459784,-0.0202816203,-0.018314434,-0.016244896,-0.0185916956,-0.0210054595,-0.0232407241,-0.0249888848,-0.025995406,-0.0261576487,-0.0255004707,-0.0241055482,-0.0221200674,-0.0197885048,-0.0171592318,-0.0198319369,-0.0226388465,-0.0252810211,-0.0273570195,-0.0285360312,-0.0287095817,-0.0279328849,-0.0262966624,-0.0239593624,-0.0212285575,-0.0179607029,-0.0209482254,-0.0241490093,-0.0272091213,-0.0296197375,-0.0309637309,-0.0311398656,-0.0302467227,-0.0283739362,-0.0256781184,-0.0225404433,-0.0185915946,-0.0218492978,-0.0253995346,-0.0288369998,-0.03154319,-0.0330199631,-0.0331911396,-0.0322027634,-0.030126953,-0.027106884,-0.0236030828,-0.0189964472,-0.0224408559,-0.02623859,-0.0299425935,-0.0328452722,-0.0343976892,-0.0345628703,-0.0335193906,-0.0313056202,-0.0280511489,-0.0242875923,-0.0191333693,-0.022646958,-0.0265362038,-0.0303293881,-0.0332793763,-0.0348373993,-0.0350025804,-0.0339534946,-0.0316924148,-0.0283487627,-0.0244936944,-0.0189851722,-0.0224342553,-0.0262327531,-0.0299096909,-0.0327467946,-0.0342453379,-0.0344165144,-0.0334063679,-0.0311996441,-0.0279401024,-0.0241880403,-0.0185651856,-0.0218262874,-0.0253699746,-0.0287553899,-0.0313533078,-0.0327428873,-0.0329190221,-0.031980293,-0.0299202048,-0.0268990837,-0.0234185053,-0.0179141714,-0.0208969536,-0.0240765425,-0.0270660171,-0.029354498,-0.0306026231,-0.0307761736,-0.0299303634,-0.0280816583,-0.0253970584,-0.0222935742,-0.0170895969,-0.019746218,-0.0225167557,-0.02507942,-0.0270397245,-0.0281300003,-0.028292243,-0.0275513104,-0.0259442441,-0.0236313637,-0.0209430272,-0.016244896,-0.0185916956,-0.0210054595,-0.0232407241,-0.0249888848,-0.025995406,-0.0261576487,-0.0255004707,-0.0241055482,-0.0221200674,-0.0197885048,-0.0173799094,-0.0201355046,-0.0230453051,-0.0257977267,-0.0279651816,-0.0291921075,-0.0293673946,-0.0285580184,-0.0268588938,-0.024428988,-0.0215920922,-0.0184618564,-0.0216577146,-0.0251287088,-0.0284901164,-0.0311537248,-0.0326219703,-0.0327931468,-0.0318132982,-0.0297800696,-0.0268360582,-0.0234114996,-0.0194265505,-0.0230635171,-0.0271266755,-0.0311584286,-0.0343674632,-0.0360715111,-0.0362130939,-0.0350628536,-0.0326960968,-0.0291931324,-0.0251272646,-0.0201985929,-0.0242286888,-0.0288473939,-0.0335313061,-0.0372582752,-0.0391473454,-0.0392389169,-0.0379534023,-0.0352979515,-0.0312493882,-0.0265638651,-0.0207018568,-0.025014513,-0.0300515412,-0.0352324715,-0.0393195609,-0.0412941837,-0.0413438743,-0.0399989692,-0.0371438384,-0.0326654918,-0.0275125722,-0.0208762135,-0.0253004136,-0.0305084578,-0.0358734197,-0.0400371656,-0.0419840826,-0.0420337732,-0.0407165739,-0.0377847866,-0.0331224084,-0.0277984728,-0.0206956885,-0.0250281992,-0.0300951261,-0.0352495035,-0.0391834275,-0.0410233733,-0.0411149448,-0.0398785546,-0.0370161489,-0.0324971204,-0.0273633755,-0.0201787715,-0.0242327286,-0.028878958,-0.0334988852,-0.0369922959,-0.0386835377,-0.0388251205,-0.0376876862,-0.0350365534,-0.0309454148,-0.0262964761,-0.0193835683,-0.0230311578,-0.027092018,-0.0310248798,-0.0339967034,-0.0355096325,-0.0356808089,-0.0346562768,-0.032314833,-0.0287993674,-0.0247849428,-0.0183900298,-0.0215759907,-0.0250117043,-0.0282576378,-0.0307194544,-0.0320323573,-0.0322076445,-0.0313122912,-0.0293188049,-0.0263953872,-0.0230325784,-0.0171592318,-0.0198319369,-0.0226388465,-0.0252810211,-0.0273570195,-0.0285360312,-0.0287095817,-0.0279328849,-0.0262966624,-0.0239593624,-0.0212285575,-0.0184618564,-0.0216577146,-0.0251287088,-0.0284901164,-0.0311537248,-0.0326219703,-0.0327931468,-0.0318132982,-0.0297800696,-0.0268360582,-0.0234114996,-0.0197252269,-0.0235055223,-0.0277680251,-0.0320399715,-0.0354596958,-0.0372556567,-0.037376673,-0.0361588426,-0.0336836982,-0.0299897372,-0.0256969861,-0.0208722052,-0.0252618974,-0.0304097986,-0.0357749164,-0.0401201284,-0.0422428634,-0.0422329755,-0.0407660577,-0.0378599966,-0.0332765749,-0.0279460699,-0.021806702,-0.0267633333,-0.0327988478,-0.0393416087,-0.044677489,-0.0470356788,-0.0468171622,-0.045158178,-0.0419109435,-0.0363522754,-0.0299122087,-0.0224262695,-0.0278098013,-0.0345713421};
  
  double pot2DZ[121] = {-0.0235165772,-0.025995406,-0.0285360312,-0.0309637309,-0.0330199631,-0.0343976892,-0.0348373993,-0.0342453379,-0.0327428873,-0.0306026231,-0.0281300003,-0.025995406,-0.0291921075,-0.0326219703,-0.0360715111,-0.0391473454,-0.0412941837,-0.0419840826,-0.0410233733,-0.0386835377,-0.0355096325,-0.0320323573,-0.0285360312,-0.0326219703,-0.0372556567,-0.0422428634,-0.0470356788,-0.050597752,-0.0517356294,-0.0500156815,-0.0461475027,-0.0412961,-0.0363639461,-0.0309637309,-0.0360715111,-0.0422428634,-0.0494927772,-0.0572916284,-0.0637454481,-0.0657190929,-0.0621567717,-0.0553274598,-0.0477875008,-0.0408636339,-0.0330199631,-0.0391473454,-0.0470356788,-0.0572916284,-0.0702972812,-0.0835042993,-0.0870245225,-0.0781417261,-0.0656256899,-0.0542276016,-0.0449563637,-0.0343976892,-0.0412941837,-0.050597752,-0.0637454481,-0.0835042993,0,0,-0.0945330902,-0.0741852171,-0.0589697212,-0.047760186,-0.0348373993,-0.0419840826,-0.0517356294,-0.0657190929,-0.0870245225,0,0,-0.0980533134,-0.0761588619,-0.0601075987,-0.0484500849,-0.0342453379,-0.0410233733,-0.0500156815,-0.0621567717,-0.0781417261,-0.0945330902,-0.0980533134,-0.085986171,-0.0704908333,-0.0572076043,-0.0468323916,-0.0327428873,-0.0386835377,-0.0461475027,-0.0553274598,-0.0656256899,-0.0741852171,-0.0761588619,-0.0704908333,-0.0611621424,-0.0516921401,-0.0434756605,-0.0306026231,-0.0355096325,-0.0412961,-0.0477875008,-0.0542276016,-0.0589697212,-0.0601075987,-0.0572076043,-0.0516921401,-0.0453365432,-0.0392516082,-0.0281300003,-0.0320323573,-0.0363639461,-0.0408636339,-0.0449563637,-0.047760186,-0.0484500849,-0.0468323916,-0.0434756605,-0.0392516082,-0.0348726072};
  
  string dx[9] = {"# Data from PBAM Electrostat run",
    "# My runname is "+test_dir_loc+"test.dx and "
    + "units internal",
    "object 1 class gridpositions counts 11 11 11",
    "origin -7.2 -7.2 -15.2",
    "delta 1.30909 0.0e+00 0.0e+00",
    "delta 0.0e00 1.30909 0.0e+00",
    "delta 0.0e00 0.0e+00 2.03636",
    "object 2 class gridconnections counts 11 11 11",
    "object 3 class array type double rank 0 items 1331 data follows"};
  
  string potx[7] = {"# Data from PBAM Electrostat run",
    "# My runname is "+test_dir_loc+"pot_x_1.00.dat",
    "units internal", "grid 11 11", "axis x 0.654545",
    "origin -7.2 -15.2", "delta 1.30909 2.03636"};
  string poty[7] = {"# Data from PBAM Electrostat run",
    "# My runname is "+test_dir_loc+"pot_y_4.00.dat",
    "units internal", "grid 11 11", "axis y 4.58182",
    "origin -7.2 -15.2", "delta 1.30909 2.03636"};
  string potz[7] = {"# Data from PBAM Electrostat run",
    "# My runname is "+test_dir_loc+"pot_z_-1.00.dat",
    "units internal","grid 11 11","axis z -0.945455",
    "origin -7.2 -7.2","delta 1.30909 1.30909"};
  
  
} ; // end EnergyForceUTest

TEST_F(ElecUTest, detailsCheck)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  vector<double>  mins = EstatTest.get_mins();
  vector<double>  maxs = EstatTest.get_maxs();
  vector<int>     npts = EstatTest.get_npts();
  vector<double>  bins = EstatTest.get_bins();
  
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR( mins[i]/min[i], 1, preclim);
    EXPECT_NEAR(    maxs[i]/7.2, 1, preclim);
    EXPECT_NEAR(   npts[i]/11.0, 1, preclim);
    EXPECT_NEAR( bins[i]/bin[i], 1, preclim);
  }
}

TEST_F(ElecUTest, detailsCheck2)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 20);
  vector<double>  mins = EstatTest.get_mins();
  vector<double>  maxs = EstatTest.get_maxs();
  vector<int>     npts = EstatTest.get_npts();
  vector<double>  bins = EstatTest.get_bins();
  
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR( mins[i]/min[i], 1, preclim);
    EXPECT_NEAR(    maxs[i]/7.2, 1, preclim);
    EXPECT_NEAR(   npts[i]/20.0, 1, preclim);
    EXPECT_NEAR( bins[i]/bn2[i], 1, preclim);
  }
}

TEST_F(ElecUTest, checkPot)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  vector<vector<vector<double > > > esp = EstatTest.get_potential();
  
  int ct = 0;
  for (int i = 0; i < esp.size(); i++)
    for (int j = 0; j < esp[0].size(); j++)
      for (int k = 0; k < esp[0][0].size(); k++)
      {
        if (ct < 300)  EXPECT_NEAR( esp[i][j][k]/potE[ct], 1, preclim);
        else           return;
        ct++;
      }
}

TEST_F(ElecUTest, printDX)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  EstatTest.print_dx(test_dir_loc+"test.dx");
  
  string inputLine;
  ifstream fin(test_dir_loc+"test.dx");
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 9) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == dx[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, checkPOTZ)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = -5.4;
  char pot[50];
  sprintf(pot, "pot_z_%.2f.dat", val);
  EstatTest.print_grid("z", val, test_dir_loc+string(pot));
  
  vector<vector<double > > p2d = EstatTest.get_pot2d();
  
  int ct = 0;
  for (int i = 0; i < p2d.size(); i++)
    for (int j = 0; j < p2d[0].size(); j++)
    {
      if ((ct < 121) && (pot2DZ[ct] != 0))
        EXPECT_NEAR( p2d[i][j]/pot2DZ[ct], 1, preclim);
      else
        return;
      ct++;
    }
}


TEST_F(ElecUTest, checkGridOutRange)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-4, 1000); ASolvTest->solve_gradA(1E-4, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = 100.8;
  char pot[50];
  sprintf(pot, "pot_x_%.2f.dat", val);
  try
  {
    EstatTest.print_grid("x", val, test_dir_loc+string(pot));
    FAIL();
  }
  catch( const ValueOutOfRange& err )
  {
    // check exception
    string exp = "x value 100.800000 out of range. It is greater than 7.200000";
    EXPECT_EQ(string(err.what()), exp);
  }
}

TEST_F(ElecUTest, checkGridOutRange2)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-4, 1000); ASolvTest->solve_gradA(1E-4, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = -50.4;
  char pot[50];
  sprintf(pot, "pot_z_%.2f.dat", val);
  try
  {
    EstatTest.print_grid("z", val, test_dir_loc+string(pot));
    FAIL();
  }
  catch( const ValueOutOfRange& err )
  {
    // check exception
    string exp = "z value -50.400000 out of range. It is less than -15.200000";
    EXPECT_EQ(string(err.what()), exp);
  }
}

TEST_F(ElecUTest, printPOTX)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = 1;
  char pot[50];
  sprintf(pot, "pot_x_%.2f.dat", val);
  EstatTest.print_grid("x", val, test_dir_loc+string(pot));
  
  string inputLine;
  ifstream fin(test_dir_loc+pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == potx[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOTY)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = 4;
  char pot[50];
  sprintf(pot, "pot_y_%.2f.dat", val);
  EstatTest.print_grid("y", val, test_dir_loc+string(pot));
  
  string inputLine;
  ifstream fin(test_dir_loc+pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == poty[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOTZ)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12, 1000); ASolvTest->solve_gradA(1E-12, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = -1;
  char pot[50];
  sprintf(pot, "pot_z_%.2f.dat", val);
  EstatTest.print_grid("z", val, test_dir_loc+string(pot));
  
  string inputLine;
  ifstream fin(test_dir_loc+pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == potz[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOT)
{
  vector<Molecule> mol_sing_;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
            Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
            Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < 9; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol_sing_.push_back( molNew );
  }
  
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_sing_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-3, 1000); ASolvTest->solve_gradA(1E-3, 1000);
  
  Electrostatic EstatTest( ASolvTest, 11);
  double val = 0; char pot[50];
  sprintf(pot, "pot_z_%.2f.dat", val);
  EstatTest.print_grid("z", val, test_dir_loc+string(pot));
  
  EXPECT_TRUE( 0 == 0);
}



#endif /* ElectrostaticsUnitTest_h */

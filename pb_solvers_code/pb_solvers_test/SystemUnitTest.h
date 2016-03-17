//
//  SystemUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/9/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef SystemUnitTest_h
#define SystemUnitTest_h

#include "System.h"

class MoleculeUTest : public ::testing::Test
{
public :
  
protected :
  virtual void SetUp() {}
  virtual void TearDown() {}
};


TEST_F(MoleculeUTest, checkUserSpecRadCent)
{
  Pt pos(0.0,0.0,-5.0);
  int M = 3; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=0; posCharges[0] = pos;
  charges[1]=2.0; vdW[1]=0; posCharges[1] = pos + Pt(1.0, 0.0, 0.0);
  charges[2]=2.0; vdW[2]=0; posCharges[2] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos);
  
  ASSERT_EQ(      3, molNew.get_m());
  ASSERT_EQ(    2.0, molNew.get_a());
  ASSERT_EQ( "stat", molNew.get_type());
  ASSERT_EQ(    0.0, molNew.get_drot());
  ASSERT_EQ(    0.0, molNew.get_dtrans());
  
  ASSERT_EQ( 0.0, molNew.get_center().x());
  ASSERT_EQ( 0.0, molNew.get_center().y());
  ASSERT_EQ(-5.0, molNew.get_center().z());
  
  ASSERT_EQ( 0.0, molNew.get_posj(0).x());
  ASSERT_EQ( 0.0, molNew.get_posj(0).y());
  
  ASSERT_EQ( 1.0, molNew.get_posj(1).x());
  ASSERT_EQ( 0.0, molNew.get_posj(1).z());
  
  ASSERT_EQ( 0.0, molNew.get_posj(2).x());
  ASSERT_EQ( 1.0, molNew.get_posj(2).y());
  
  ASSERT_EQ( 0.0, molNew.get_posj_realspace(0).x());
  ASSERT_EQ( 0.0, molNew.get_posj_realspace(0).y());
  ASSERT_EQ(-5.0, molNew.get_posj_realspace(0).z());
  
  ASSERT_EQ( 1.0, molNew.get_posj_realspace(1).x());
  ASSERT_EQ( 0.0, molNew.get_posj_realspace(1).y());
  ASSERT_EQ(-5.0, molNew.get_posj_realspace(1).z());
}


TEST_F(MoleculeUTest, checkUserSpecRad)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "move", 13.7, charges, posCharges, vdW, 0.04, 0.34);
  
  ASSERT_EQ(      2, molNew.get_m());
  ASSERT_EQ(   13.7, molNew.get_a());
  ASSERT_EQ( "move", molNew.get_type());
  ASSERT_EQ(   0.04, molNew.get_drot());
  ASSERT_EQ(   0.34, molNew.get_dtrans());
  
  ASSERT_EQ( -9.5, molNew.get_center().x());
  ASSERT_EQ( 23.9, molNew.get_center().y());
  ASSERT_EQ( -8.7, molNew.get_center().z());

  ASSERT_EQ( 0.5, molNew.get_posj(0).x());
  ASSERT_EQ(-0.5, molNew.get_posj(0).y());
  ASSERT_EQ( 0.0, molNew.get_posj(0).z());
  
  ASSERT_EQ(-0.5, molNew.get_posj(1).x());
  ASSERT_EQ( 0.5, molNew.get_posj(1).y());
  ASSERT_EQ( 0.0, molNew.get_posj(1).z());
  
  ASSERT_EQ( -9.0, molNew.get_posj_realspace(0).x());
  ASSERT_EQ( 23.4, molNew.get_posj_realspace(0).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(0).z());
  
  ASSERT_EQ(-10.0, molNew.get_posj_realspace(1).x());
  ASSERT_EQ( 24.4, molNew.get_posj_realspace(1).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(1).z());
}

TEST_F(MoleculeUTest, checkUserSpecCent)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, pos, 0.24);
  
  ASSERT_EQ(    2, molNew.get_m());
  ASSERT_EQ("rot", molNew.get_type());
  ASSERT_EQ( 0.24, molNew.get_drot());
  ASSERT_EQ( 0.00, molNew.get_dtrans());
  EXPECT_NEAR( 7.32, molNew.get_a(), preclim);
  
  ASSERT_EQ(-10.0, molNew.get_center().x());
  ASSERT_EQ( 23.4, molNew.get_center().y());
  ASSERT_EQ( -8.7, molNew.get_center().z());
  
  ASSERT_EQ( 1.0, molNew.get_posj(0).x());
  ASSERT_EQ( 0.0, molNew.get_posj(0).y());
  ASSERT_EQ( 0.0, molNew.get_posj(0).z());
  
  ASSERT_EQ( 0.0, molNew.get_posj(1).x());
  ASSERT_EQ( 1.0, molNew.get_posj(1).y());
  ASSERT_EQ( 0.0, molNew.get_posj(1).z());
  
  ASSERT_EQ( -9.0, molNew.get_posj_realspace(0).x());
  ASSERT_EQ( 23.4, molNew.get_posj_realspace(0).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(0).z());
  
  ASSERT_EQ(-10.0, molNew.get_posj_realspace(1).x());
  ASSERT_EQ( 24.4, molNew.get_posj_realspace(1).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(1).z());
}

TEST_F(MoleculeUTest, checkCreateCen)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, 0.24);
  
  ASSERT_EQ(    2, molNew.get_m());
  ASSERT_EQ("rot", molNew.get_type());
  ASSERT_EQ( 0.24, molNew.get_drot());
  ASSERT_EQ( 0.00, molNew.get_dtrans());
  EXPECT_NEAR( 7.0271067811865, molNew.get_a(), preclim);
  
  ASSERT_EQ( -9.5, molNew.get_center().x());
  ASSERT_EQ( 23.9, molNew.get_center().y());
  ASSERT_EQ( -8.7, molNew.get_center().z());
  
  ASSERT_EQ( 0.5, molNew.get_posj(0).x());
  ASSERT_EQ(-0.5, molNew.get_posj(0).y());
  ASSERT_EQ( 0.0, molNew.get_posj(0).z());
  
  ASSERT_EQ(-0.5, molNew.get_posj(1).x());
  ASSERT_EQ( 0.5, molNew.get_posj(1).y());
  ASSERT_EQ( 0.0, molNew.get_posj(1).z());
  
  ASSERT_EQ( -9.0, molNew.get_posj_realspace(0).x());
  ASSERT_EQ( 23.4, molNew.get_posj_realspace(0).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(0).z());
  
  ASSERT_EQ(-10.0, molNew.get_posj_realspace(1).x());
  ASSERT_EQ( 24.4, molNew.get_posj_realspace(1).y());
  ASSERT_EQ( -8.7, molNew.get_posj_realspace(1).z());
}

TEST_F(MoleculeUTest, translate)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, 0.24);
  molNew.translate( Pt( 3.0, -4.5, 10.21));
  
  EXPECT_NEAR( -6.5, molNew.get_center().x(), preclim);
  EXPECT_NEAR( 19.4, molNew.get_center().y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_center().z(), preclim);
  
  EXPECT_NEAR( 0.5, molNew.get_posj(0).x(), preclim);
  EXPECT_NEAR(-0.5, molNew.get_posj(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj(0).z(), preclim);
  
  EXPECT_NEAR(-0.5, molNew.get_posj(1).x(), preclim);
  EXPECT_NEAR( 0.5, molNew.get_posj(1).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj(1).z(), preclim);
  
  EXPECT_NEAR( -6.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 18.9, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -7.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( 19.9, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_posj_realspace(1).z(), preclim);
}

TEST_F(MoleculeUTest, rotateSimple)
{
  Pt pos( 0.0, 0.0, 0.0);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, pos, 0.24);
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  ASSERT_EQ( 0.0, molNew.get_center().x());
  ASSERT_EQ( 0.0, molNew.get_center().y());
  ASSERT_EQ( 0.0, molNew.get_center().z());
  
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 1.0, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
  
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  EXPECT_NEAR(-1.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
}

TEST_F(MoleculeUTest, rotate2)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, pos, 0.24);
  molNew.rotate( Quat( 1.0, Pt(1.0, 1.0, 1.0)));
  
  ASSERT_EQ(-10.0, molNew.get_center().x());
  ASSERT_EQ( 23.4, molNew.get_center().y());
  ASSERT_EQ( -8.7, molNew.get_center().z());
  
  EXPECT_NEAR( 0.693534871, molNew.get_posj(0).x(), preclim);
  EXPECT_NEAR( 0.639056064, molNew.get_posj(0).y(), preclim);
  EXPECT_NEAR(-0.332590935, molNew.get_posj(0).z(), preclim);
  
  EXPECT_NEAR(-0.332590935, molNew.get_posj(1).x(), preclim);
  EXPECT_NEAR( 0.693534871, molNew.get_posj(1).y(), preclim);
  EXPECT_NEAR( 0.639056064, molNew.get_posj(1).z(), preclim);
  
  EXPECT_NEAR( -9.30646513, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 24.03905610, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( -9.03259093, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -10.3325909, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( 24.09353490, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR( -8.06094394, molNew.get_posj_realspace(1).z(), preclim);
}

class SystemUTest : public ::testing::Test
{
public :
  
protected :
  Constants const_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
};


TEST_F(SystemUTest, checkOverlap)
{
  vector < Molecule > mol_;
  Pt pos[2] = { Pt( 0.0, 0.0, -5.0), Pt( 0.0, 0.0, 0.0)};
  double rad[2] = { 5.0, 3.7 };
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  try
  {
    System sys( mol_ );
    FAIL();
  }
  catch( const OverlappingMoleculeException& err )
  {
    // check exception
    string error_exp = "Molecule 0 & 1 overlap";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(SystemUTest, checkPBCOverlap)
{
  vector < Molecule > mol_; const int nMol = 3;
  Pt pos[nMol] = { Pt( 0.0, 0.0, -5.0), Pt(10.0,7.8,25.0), Pt(-10.0,7.8,25.0) };
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0.1; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0.1; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0.1; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  try
  {
    System sys( mol_, 20.0, 10.0 );
    FAIL();
  }
  catch( const OverlappingMoleculeException& err )
  {
    // check exception
    string error_exp = "Molecule 1 & 2 overlap";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(SystemUTest, checkVals)
{
  vector < Molecule > mol_;
  double cutoff = 45.876;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,7.8,25.0) };
  double rad[3] = { 5.0, 3.7, 8.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  System sys( mol_, cutoff );
  EXPECT_NEAR( 5.7666666667, sys.get_lambda(), preclim);
  EXPECT_NEAR( cutoff/sys.get_cutoff(), 1.0, preclim);
  
  sys.set_time( 1.435);
  EXPECT_NEAR( 1.435/sys.get_time(), 1.0, preclim);
  
  EXPECT_EQ(sys.less_than_cutoff(Pt( 5.3, 0.34,  -2.13)), true);
  EXPECT_EQ(sys.less_than_cutoff(Pt(25.3,79.34, -12.13)), false);

  sys.translate_mol(0, Pt(0.0, 0.0, -5.0));
  
  EXPECT_NEAR(  0.0, sys.get_centeri(0).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(0).y(), preclim);
  EXPECT_NEAR(-10.0, sys.get_centeri(0).z(), preclim);
  
  sys.translate_mol(1, Pt(-10.0, -7.8, -25.0));
  
  EXPECT_NEAR(  0.0, sys.get_centeri(1).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(1).y(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(1).z(), preclim);
  
  sys.rotate_mol(1, Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  ASSERT_EQ( 0.0, sys.get_centeri(1).x());
  ASSERT_EQ( 0.0, sys.get_centeri(1).y());
  ASSERT_EQ( 0.0, sys.get_centeri(1).z());
  
  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).x(), preclim);
  EXPECT_NEAR( 1.0, sys.get_posij(1, 1).y(), preclim);
  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).z(), preclim);
  
  EXPECT_NEAR( -1.0, sys.get_posij(1, 2).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).y(), preclim);
  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).z(), preclim);
}

TEST_F(SystemUTest, changeCutoff)
{
  vector < Molecule > mol_;
  double cutoff = 45.876;
  double boxl   = 55.876;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,17.8,15.0) };
  double rad[3] = { 5.0, 3.7, 8.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  System sys( mol_, cutoff, boxl );
  EXPECT_NEAR( (boxl/2.0)/sys.get_cutoff(), 1.0, preclim);
}

TEST_F(SystemUTest, makeSubsystem)
{
  vector <int> subInd(2);
  vector < Molecule > mol_;
  double cutoff =  45.876;
  double boxl   = 105.876;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,7.8,25.0) };
  double rad[3] = { 5.0, 3.7, 8.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "rot", rad[molInd], chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  subInd[0] = 1; subInd[1] = 2;
  System sys( mol_, cutoff, boxl );
  sys.set_time( 1.435);
  
  System sub = sys.get_subsystem(subInd);
  
  EXPECT_NEAR( 6.15, sub.get_lambda(), preclim);
  EXPECT_NEAR( cutoff/sub.get_cutoff(), 1.0, preclim);
  
  EXPECT_NEAR( 1.435/sub.get_time(), 1.0, preclim);
  
  EXPECT_EQ(sub.less_than_cutoff(Pt( 5.3, 0.34,  -2.13)), true);
  EXPECT_EQ(sub.less_than_cutoff(Pt(25.3,79.34, -12.13)), false);
  
  EXPECT_EQ(sub.get_typei(0), "rot");
  EXPECT_EQ(sub.get_n(), 2);
  
}

TEST_F(SystemUTest, PBCcheck)
{
  vector < Molecule > mol_;
  double cutoff =  5.00;
  double boxl   = 20.00;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,17.8,15.0) };
  double rad[3] = { 1.0, 1.7, 1.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "rot", rad[molInd], chg, poschg, vdW, pos[molInd]);
    mol_.push_back( molNew );
  }

  System sys( mol_, cutoff, boxl );
  Pt dis01 = sys.get_pbc_dist_vec(0, 1);
  EXPECT_NEAR(  10/dis01.x(), 1.0, preclim);
  EXPECT_NEAR(-7.8/dis01.y(), 1.0, preclim);
  EXPECT_NEAR(  10/dis01.z(), 1.0, preclim);
  
  Pt dis02 = sys.get_pbc_dist_vec(0, 2);
  EXPECT_NEAR( -10/dis02.x(), 1.0, preclim);
  EXPECT_NEAR( 2.2/dis02.y(), 1.0, preclim);
  EXPECT_NEAR(   0,     dis02.z(), preclim);
  
  Pt dis12 = sys.get_pbc_dist_vec(1, 2);
  EXPECT_NEAR( 0.0,  dis12.x(), preclim);
  EXPECT_NEAR(  10/dis12.y(), 1.0, preclim);
  EXPECT_NEAR( -10/dis12.z(), 1.0, preclim);

}



#endif /* SystemUnitTest_h */

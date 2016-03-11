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

class SystemUTest : public ::testing::Test
{
public :
  
protected :
  virtual void SetUp() {}
  virtual void TearDown() {}
};


TEST_F(SystemUTest, checkVals)
{
  
}


#endif /* SystemUnitTest_h */

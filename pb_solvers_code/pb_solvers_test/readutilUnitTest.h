//
//  readutilUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef readutilUnitTest_h
#define readutilUnitTest_h

#include "readutil.h"

/*
 Class to test opening PQR and XYZ files
 */
class ReadUtilUTest : public ::testing::Test
{
public :
  
protected :
  // This will change according to where you're running tests!
  string test_dir_ = "/Users/lfelberg/PBSAM/pb_solvers/pb_solvers_code/test/";
  
};

TEST_F(ReadUtilUTest, checkPQRExceptions)
{
  string path = test_dir_ + "none.pqr";
  try
  {
    PQRFile PQRtest(path, 10);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readPQR)
{
  string PQR = test_dir_ + "2charge.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  vector<Pt> my_cents = PQRtest.get_cg_centers();
  
  ASSERT_EQ(true, PQRtest.get_cg());
  ASSERT_EQ(   2, PQRtest.get_M());
  ASSERT_EQ( PQR, PQRtest.get_path());
  ASSERT_EQ( my_atoms[0].x(), 0);
  ASSERT_EQ( my_atoms[0].y(), 0);
  ASSERT_EQ( my_atoms[0].z(), 0);
  
  ASSERT_EQ( my_atoms[1].x(), 0);
  ASSERT_EQ( my_atoms[1].y(), 1);
  ASSERT_EQ( my_atoms[1].z(), 0);
  
  ASSERT_EQ(   -1, PQRtest.get_charges()[0]);
  ASSERT_EQ( 1.67, PQRtest.get_charges()[1]);
  ASSERT_EQ( 0.45, PQRtest.get_radii()[0]);
  ASSERT_EQ( 0.20, PQRtest.get_radii()[1]);
  
  ASSERT_EQ(    1, PQRtest.get_cg_centers().size());
  ASSERT_EQ(    0, my_cents[0].x());
  ASSERT_EQ(    0, my_cents[0].y());
  ASSERT_EQ(    0, my_cents[0].z());
  ASSERT_EQ( 1.85, PQRtest.get_cg_radii()[0]);
}

TEST_F(ReadUtilUTest, readPQRNoCen)
{
  string PQR = test_dir_ + "2charge_nocen.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  
  ASSERT_EQ(false, PQRtest.get_cg());
  ASSERT_EQ(   4, PQRtest.get_M());
  ASSERT_EQ( PQR, PQRtest.get_path());
  ASSERT_EQ( my_atoms[0].x(), 0);
  ASSERT_EQ( my_atoms[0].y(), 0);
  ASSERT_EQ( my_atoms[0].z(), 0);
  
  ASSERT_EQ( my_atoms[2].x(), 5);
  ASSERT_EQ( my_atoms[2].y(), 0);
  ASSERT_EQ( my_atoms[2].z(), 0);
  
  ASSERT_EQ(   -1, PQRtest.get_charges()[0]);
  ASSERT_EQ( 7.56, PQRtest.get_charges()[1]);
  ASSERT_EQ( 0.50, PQRtest.get_radii()[0]);
  ASSERT_EQ( 0.87, PQRtest.get_radii()[1]);
}


TEST_F(ReadUtilUTest, checkXYZExceptions)
{
  string path = test_dir_ + "none.xyz";
  try
  {
    PQRFile XYZtest(path, 10);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, checkXYZExceptionShort)
{
  string xyz = test_dir_ + "2charge_short.xyz";
  try
  {
    XYZFile XYZtest(xyz, 10);
    FAIL();
  }
  catch( const NotEnoughCoordsException& err )
  {
    // check exception
    string error_exp = "File has 5 lines, need 10";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}


TEST_F(ReadUtilUTest, readXYZ)
{
  string xyz = test_dir_ + "2charge.xyz";
  XYZFile XYZtest(xyz, 10);
  
  ASSERT_EQ(  10, XYZtest.get_nmols());
  ASSERT_EQ( xyz, XYZtest.get_path());
}


#endif /* readutilUnitTest_h */

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

  vector<double> sp_x = {2.965, 4.52, 7.092, 7.782, 7.614, 4.82, 4.093, 7.83, 
    6.265, 12.017, 9.993, 8.764, 12.731, 8.46, 4.095, 6.813, 8.422};
  vector<double> sp_z = {-1.084, -3.482, -6.049, 1.145, 4.029, 2.408, -0.923, 
    2.59, -6.542, 3.608, -7.069, 3.456, -4.511, -2.487, -3.01, -0.62, 1.101};
  vector<double> np_y = {-0.636, -0.442, -0.574, -0.594, -0.524, 0.526, -0.048,
   0.512, -0.779, -0.442, 0.201, -0.895, 0.184, -0.595, -0.171, -0.111, -0.75};
  vector<double> np_z = {-0.175, -0.581, -0.81, -0.453, -0.01, 0.656, -0.377,
   0.758, -0.608, 0.723, -0.842, 0.091, -0.82, -0.382, -0.296, -0.774, -0.22};
};


// MSMS section
TEST_F(ReadUtilUTest, checkMSMSExceptions)
{
  string path = test_dir_loc + "none.face";
  try
  {
    MSMSFile MSMStest(path);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readMSMS)
{
  string path = test_dir_loc + "test.vert";
  MSMSFile MSMStest(path);
  vector<Pt> my_sp = MSMStest.get_sp();
  vector<Pt> my_np = MSMStest.get_np();

  ASSERT_EQ( my_sp.size(), 1669);
  ASSERT_EQ( my_np.size(), 1669);

  int ct = 0;
  for (int i=0; i<sp_x.size(); i += 100)
  {
    EXPECT_NEAR( my_sp[i].x(), sp_x[ct], preclim);
    EXPECT_NEAR( my_sp[i].z(), sp_z[ct], preclim);
    EXPECT_NEAR( my_np[i].y(), np_y[ct], preclim);
    EXPECT_NEAR( my_np[i].z(), np_z[ct], preclim);
    ct++;
  }

}

// Contact section
TEST_F(ReadUtilUTest, checkContactExceptions)
{
  string path = test_dir_loc + "none.cont";
  try
  {
    ContactFile contTest(path);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readContact)
{
  string path = test_dir_loc + "test.cont";
  ContactFile contTest(path);
  vector<vector<int > > my_pr = contTest.get_at_pairs();

  ASSERT_EQ( contTest.get_moltype1(), 0);
  ASSERT_EQ( contTest.get_moltype2(), 21);

  vector<int > pair1 = {22, 47, 17};
  vector<int > pair2 = {46, 108, 455};
  vector<double > dst = {3.0, 2.57, 9.57};

  for (int i=0; i<my_pr.size(); i++)
  {
    ASSERT_EQ( my_pr[i][0], pair1[i]);
    ASSERT_EQ( my_pr[i][1], pair2[i]);
    EXPECT_NEAR( contTest.get_dists()[i], dst[i], preclim);
  }

}

// PQR section
TEST_F(ReadUtilUTest, checkPQRExceptions)
{
  string path = test_dir_loc + "none.pqr";
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
  string PQR = test_dir_loc + "2charge.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  vector<Pt> my_cents = PQRtest.get_cg_centers();
  
//  ASSERT_EQ(true, PQRtest.get_cg());
//  ASSERT_EQ(   2, PQRtest.get_M());  // !!TODO: Change
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
  string PQR = test_dir_loc + "2charge_nocen.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  
//  ASSERT_EQ(false, PQRtest.get_cg());
//  ASSERT_EQ(   4, PQRtest.get_M());  // !!Change
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
  string path = test_dir_loc + "none.xyz";
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
  string xyz = test_dir_loc + "2charge_short.xyz";
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
  string xyz = test_dir_loc + "2charge.xyz";
  XYZFile XYZtest(xyz, 10);
  
  ASSERT_EQ(  10, XYZtest.get_nmols());
  ASSERT_EQ( xyz, XYZtest.get_path());
}


#endif /* readutilUnitTest_h */

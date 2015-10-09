//
//  utilUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/9/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef utilUnitTest_h
#define utilUnitTest_h

#include "util.h"

/*
 Class for testing euclidean points
 */
class EuPointUTest : public ::testing::Test
{
  public :
  
  protected :
  
  EuPoint <double> testPt_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};

TEST_F(EuPointUTest, setGetConvert)
{
  EuPoint <double> testPt( 76.4, 33.4, 47.2 );
  EXPECT_NEAR( testPt.get_x(), 76.4, preclim);
  EXPECT_NEAR( testPt.get_y(), 33.4, preclim);
  EXPECT_NEAR( testPt.get_z(), 47.2, preclim);
  
  SphPoint<double> testSpPt = testPt.convert_to_spherical();
  EXPECT_NEAR( testSpPt.get_r(),     95.814195190, preclim);
  EXPECT_NEAR( testSpPt.get_theta(),  0.412135754, preclim);
  EXPECT_NEAR( testSpPt.get_phi(),    1.055698348, preclim);
}

TEST_F(EuPointUTest, setGetConvert2)
{
  EuPoint <double> testPt( 0.0, 20.0, 5.0 );
  EXPECT_NEAR( testPt.get_x(),  0.0, preclim);
  EXPECT_NEAR( testPt.get_y(), 20.0, preclim);
  EXPECT_NEAR( testPt.get_z(),  5.0, preclim);
  
  SphPoint<double> testSpPt = testPt.convert_to_spherical();
  EXPECT_NEAR( testSpPt.get_r(),     20.615528128, preclim);
  EXPECT_NEAR( testSpPt.get_theta(),  1.570796327, preclim);
  EXPECT_NEAR( testSpPt.get_phi(),    1.325817664, preclim);
}

/*
 Class for testing spherical coordinate points
 */
class SpPointUTest : public ::testing::Test
{
  public :
  
  protected :
  
  EuPoint <double> testPt_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};

TEST_F(SpPointUTest, setGetConvert)
{
  SphPoint <double> testPt( -30.3, 0.8, 9.7 );
  EXPECT_NEAR( testPt.get_r(),     -30.3, preclim);
  EXPECT_NEAR( testPt.get_theta(),   0.8, preclim);
  EXPECT_NEAR( testPt.get_phi(),     9.7, preclim);
  
  EuPoint<double> testEuPt = testPt.convert_to_euclidean();
  EXPECT_NEAR( testEuPt.get_x(),  5.73692479, preclim);
  EXPECT_NEAR( testEuPt.get_y(),  5.90695896, preclim);
  EXPECT_NEAR( testEuPt.get_z(), 29.15965586, preclim);
}

TEST_F(SpPointUTest, setGetConvert2)
{
  SphPoint <double> testPt( 5.0, 7.5, 3.5 );
  EXPECT_NEAR( testPt.get_r(),     5.0, preclim);
  EXPECT_NEAR( testPt.get_theta(), 7.5, preclim);
  EXPECT_NEAR( testPt.get_phi(),   3.5, preclim);
  
  EuPoint<double> testEuPt = testPt.convert_to_euclidean();
  EXPECT_NEAR( testEuPt.get_x(), -0.6079693, preclim);
  EXPECT_NEAR( testEuPt.get_y(), -1.6451733, preclim);
  EXPECT_NEAR( testEuPt.get_z(), -4.6822834, preclim);
}

#endif /* utilUnitTest_h */

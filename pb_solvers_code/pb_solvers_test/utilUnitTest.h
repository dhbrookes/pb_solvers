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
  
  Point <double> testPt_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};

TEST_F(EuPointUTest, setGetConvert)
{
  Point <double> testPt( 76.4, 33.4, 47.2, false );
  EXPECT_NEAR( testPt.x(), 76.4, preclim);
  EXPECT_NEAR( testPt.y(), 33.4, preclim);
  EXPECT_NEAR( testPt.z(), 47.2, preclim);
  
  EXPECT_NEAR( testPt.r(),     95.814195190, preclim);
  EXPECT_NEAR( testPt.theta(),  1.055698348, preclim);
  EXPECT_NEAR( testPt.phi(),    0.412135754, preclim);
}

TEST_F(EuPointUTest, setGetConvert2)
{
  Point <double> testPt( 0.0, 20.0, 5.0, false );
  EXPECT_NEAR( testPt.x(),  0.0, preclim);
  EXPECT_NEAR( testPt.y(), 20.0, preclim);
  EXPECT_NEAR( testPt.z(),  5.0, preclim);
  
  EXPECT_NEAR( testPt.r(),     20.615528128, preclim);
  EXPECT_NEAR( testPt.theta(),  1.325817664, preclim);
  EXPECT_NEAR( testPt.phi(),    1.570796327, preclim);
}

TEST_F(EuPointUTest, setGetConvert3)
{
  Point <double> testPt( 0.0, 0.0,-5.0, false );
  EXPECT_NEAR( testPt.x(),  0.0, preclim);
  EXPECT_NEAR( testPt.y(),  0.0, preclim);
  EXPECT_NEAR( testPt.z(), -5.0, preclim);

  EXPECT_NEAR( testPt.r(),      5.000000000, preclim);
  EXPECT_NEAR( testPt.theta(),         M_PI, preclim);
  EXPECT_NEAR( testPt.phi(),    0.000000000, preclim);
}

/*
 Class for testing spherical coordinate points
  used: http://keisan.casio.com/exec/system/1359533867 for calcs
 */
class SpPointUTest : public ::testing::Test
{
  public :
  
  protected :
  
  Point <double> testPt_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};

TEST_F(SpPointUTest, setGetConvert)
{
  Point <double> testPt( -30.3, 0.8, 9.7, true );
  EXPECT_NEAR( testPt.r(),     -30.3, preclim);
  EXPECT_NEAR( testPt.theta(),   0.8, preclim);
  EXPECT_NEAR( testPt.phi(),     9.7, preclim);
  
  EXPECT_NEAR( testPt.x(),  5.73692479, preclim);
  EXPECT_NEAR( testPt.y(),  5.90695896, preclim);
  EXPECT_NEAR( testPt.z(), 29.15965586, preclim);
}

TEST_F(SpPointUTest, setGetConvert2)
{
  Point <double> testPt( 5.0, 7.5, 3.5, true );
  EXPECT_NEAR( testPt.r(),     5.0, preclim);
  EXPECT_NEAR( testPt.theta(), 7.5, preclim);
  EXPECT_NEAR( testPt.phi(),   3.5, preclim);

  EXPECT_NEAR( testPt.x(), -0.6079693, preclim);
  EXPECT_NEAR( testPt.y(), -1.6451733, preclim);
  EXPECT_NEAR( testPt.z(), -4.6822834, preclim);
}

#endif /* utilUnitTest_h */

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
  
  EXPECT_NEAR( testPt.x(),  20.91785674, preclim);
  EXPECT_NEAR( testPt.y(),   5.90695896, preclim);
  EXPECT_NEAR( testPt.z(), -21.11021329, preclim);
}

TEST_F(SpPointUTest, setGetConvert2)
{
  Point <double> testPt( 5.0, 7.5, 3.5, true );
  EXPECT_NEAR( testPt.r(),     5.0, preclim);
  EXPECT_NEAR( testPt.theta(), 7.5, preclim);
  EXPECT_NEAR( testPt.phi(),   3.5, preclim);

  EXPECT_NEAR( testPt.x(), -4.391981755, preclim);
  EXPECT_NEAR( testPt.y(), -1.645173297, preclim);
  EXPECT_NEAR( testPt.z(),  1.733176589, preclim);
}


/*
 Class for testing euclidean points
 */
class QuatUTest : public ::testing::Test
{
  public :
  
  protected :
  
  Quaternion testQuat_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};

TEST_F(QuatUTest, InitVals) // test for constructors
{
  Quaternion testQuat;
  EXPECT_NEAR( testQuat.get_w(), 1, preclim);
  EXPECT_NEAR( testQuat.get_a(), 0, preclim);
  EXPECT_NEAR( testQuat.get_b(), 0, preclim);
  EXPECT_NEAR( testQuat.get_c(), 0, preclim);
  
  Quaternion testQuat1( 2.4, 4.6, 7.8, 9.2);
  EXPECT_NEAR( testQuat1.get_w(), 0.18278586612161321, preclim);
  EXPECT_NEAR( testQuat1.get_a(), 0.35033957673309196, preclim);
  EXPECT_NEAR( testQuat1.get_b(), 0.59405406489524293, preclim);
  EXPECT_NEAR( testQuat1.get_c(), 0.70067915346618392, preclim);
  
  Quaternion testQuat2( 2.8, Point<double>(7.6, 1.8, 11.2));
  EXPECT_NEAR( testQuat2.get_w(), 0.16996714290024101, preclim);
  EXPECT_NEAR( testQuat2.get_a(), 0.54850238458972622, preclim);
  EXPECT_NEAR( testQuat2.get_b(), 0.12990845950809307, preclim);
  EXPECT_NEAR( testQuat2.get_c(), 0.80831930360591231, preclim);
}

TEST_F(QuatUTest, Conjugate) // test for conj()
{
  Quaternion testQuat1( 2.4, 4.6, 7.8, 9.2);
  EXPECT_NEAR( testQuat1.conj().get_w(),  0.18278586612161321, preclim);
  EXPECT_NEAR( testQuat1.conj().get_a(), -0.35033957673309196, preclim);
  EXPECT_NEAR( testQuat1.conj().get_b(), -0.59405406489524293, preclim);
  EXPECT_NEAR( testQuat1.conj().get_c(), -0.70067915346618392, preclim);
}

TEST_F(QuatUTest, Multiply) // test for operator*(Quat rhs)
{
  Quaternion testQuat1( 2.4, 4.6, 7.8, 9.2);
  Quaternion testQuat2( 2.8, Point<double>(7.6, 1.8, 11.2));

  Quaternion test3 = testQuat1 * testQuat2;
  EXPECT_NEAR( test3.get_w(), -0.804639636, preclim);
  EXPECT_NEAR( test3.get_a(),  0.548965919, preclim);
  EXPECT_NEAR( test3.get_b(),  0.225853046, preclim);
  EXPECT_NEAR( test3.get_c(), -0.0134862186, preclim);
  
  test3 = testQuat2 * testQuat1;
  EXPECT_NEAR( test3.get_w(), -0.804639636, preclim);
  EXPECT_NEAR( test3.get_a(), -0.229356518, preclim);
  EXPECT_NEAR( test3.get_b(),  0.0235771586, preclim);
  EXPECT_NEAR( test3.get_c(),  0.547169774, preclim);
}

TEST_F(QuatUTest, Rotate) // test for rotate_point
{
  Quaternion testQuat1( 2.4, 4.6, 7.8, 9.2);
  Quaternion testQuat2( 2.8, Point<double>(7.6, 1.8, 11.2));
  
  Point<double> point1( 0.5, 5.0, 100.4);
  
  Point<double> test3 = testQuat1.rotate_point( point1 );
  cout << " This is my point " << test3.x() << " " << test3.y() << " " << test3.z() << endl;
  EXPECT_NEAR( test3.x()/71.5519258, 1, preclim);
  EXPECT_NEAR( test3.y()/69.9219026, 1, preclim);
  EXPECT_NEAR( test3.z()/9.83155452, 1, preclim);
  
  test3 = testQuat2.rotate_point( point1 );
    cout << " This is my point " << test3.x() << " " << test3.y() << " " << test3.z() << endl;
  EXPECT_NEAR( test3.x()/92.6298207,  1, preclim);
  EXPECT_NEAR( test3.y()/-1.96825623, 1, preclim);
  EXPECT_NEAR( test3.z()/39.0032343,  1, preclim);
}

#endif /* utilUnitTest_h */

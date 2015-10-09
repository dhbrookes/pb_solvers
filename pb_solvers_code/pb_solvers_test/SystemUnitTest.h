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


TEST_F(MoleculeUTest, checkVals)
{
  
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

/*
 
 // testing matrix multiplication
 TEST_F(MyMatrixUTest, matMul)
 {
 MyMatrix <double> testMat1( testInp_ );
 MyMatrix <double> testMat2( testInp2_ );
 MyMatrix <double> testMat3;
 
 testMat3 = testMat1 * testMat2;
 ASSERT_EQ( testMat3.get_nrows() , 2 );
 ASSERT_EQ( testMat3.get_ncols() , 2 );
 
 EXPECT_NEAR( testMat3(0,0),  50.0, preclim);
 EXPECT_NEAR( testMat3(0,1),  60.0, preclim);
 EXPECT_NEAR( testMat3(1,0), 114.0, preclim);
 EXPECT_NEAR( testMat3(1,1), 140.0, preclim);
 }
 
 // testing matrix assertions for accessing matrix members
 TEST_F(MyMatrixUTest, exceptAccess)
 {
 MyMatrix <double> testMat_( testInp_ );
 ASSERT_THROW(testMat_( 0, 5), MatrixAccessException);
 ASSERT_THROW(testMat_( 5, 0), MatrixAccessException);
 ASSERT_THROW(testMat_(-1, 0), MatrixAccessException);
 ASSERT_THROW(testMat_( 0,-1), MatrixAccessException);
 }
 
 
 */

#endif /* SystemUnitTest_h */

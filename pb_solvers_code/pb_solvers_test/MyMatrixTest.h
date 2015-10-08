//
//  MyMatrixTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/6/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef MyMatrixTest_h
#define MyMatrixTest_h

#include "MyMatrix.h"
using namespace std;

/*
 Class for testing matrix class
 */
class MyMatrixUTest  : public ::testing::Test
{
public:


protected :
  virtual void SetUp()
  {
    const int nCol = 4;
    const int nRow = 2;
    
    double val = 1.0;
    testInp_.resize(nRow);
    
    for ( int i = 0.0; i < nRow; i++ )
    {
      testInp_[i].resize(nCol);
      for ( int j = 0.0; j < nCol; j++ )
      {
        testInp_[i][j] = val;
        val += 1.0;
      }
    }
    
    testInp2_.resize(nCol);
    val = 1.0;
    
    for ( int i = 0; i < nCol; i++ )
    {
      testInp2_[i].resize(nRow);
      for ( int j = 0; j < nRow; j++ )
      {
        testInp2_[i][j] = val;
        val += 1.0;
      }
    }
    
  }
  virtual void TearDown() {}
  
  MyMatrix <double> testMat_;
  vector< vector<double> > testInp_;
  vector< vector<double> > testInp2_;
}; // end MyMatrixTest



TEST_F(MyMatrixUTest, constFromOther)
{
  MyMatrix <double> testMat_( testInp_ );
  EXPECT_NEAR( testMat_(0,0), 1.0, preclim);
  EXPECT_NEAR( testMat_(0,3), 4.0, preclim);
  EXPECT_NEAR( testMat_(1,0), 5.0, preclim);
  EXPECT_NEAR( testMat_(1,2), 7.0, preclim);
}

TEST_F(MyMatrixUTest, defaultConstruct)
{
  MyMatrix <double> testMat_( 4, 3);
  ASSERT_EQ( testMat_.get_nrows() , 4 );
  ASSERT_EQ( testMat_.get_ncols() , 3 );
  testMat_.set_val( 3, 1, 5.2);
  EXPECT_NEAR( testMat_( 3, 1), 5.2, preclim);
}

// testing matrix addition
TEST_F(MyMatrixUTest, matAdd)
{
  MyMatrix <double> testMat1( testInp_ );
  MyMatrix <double> testMat2( testInp_ );
  MyMatrix <double> testMat3;
  
  testMat3 = testMat1 + testMat2;
  
  EXPECT_NEAR( testMat3(0,0),  2.0, preclim);
  EXPECT_NEAR( testMat3(0,3),  8.0, preclim);
  EXPECT_NEAR( testMat3(1,0), 10.0, preclim);
  EXPECT_NEAR( testMat3(1,2), 14.0, preclim);
}

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



#endif /* MyMatrixTest_h */

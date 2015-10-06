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
protected:
  MyMatrix <double> testMat_;
  
public:
  MyMatrixUTest( )
  {
    const int nCol = 4;
    const int nRow = 2;
    
    double val = 1.0;
    vector< vector<double> > testInp;
    testInp.resize(nRow);
    
    for ( int i = 0.0; i < nRow; i++ )
    {
      testInp[i].resize(nCol);
      for ( int j = 0.0; j < nCol; j++ )
      {
        testInp[i][j] = val;
        val += 1.0;
      }
    }
    
    MyMatrix <double> testMat_( testInp );
    
  }

protected :
  virtual void SetUp() {}
  virtual void TearDown() {}
}; // end MyMatrixTest



TEST_F(MyMatrixUTest, settingAndCalc)
{
//  EXPECT_NEAR( ConstUTest.convert_int_to_kcal_mol( 1.0 ), 332.061203, preclim);
}




#endif /* MyMatrixTest_h */

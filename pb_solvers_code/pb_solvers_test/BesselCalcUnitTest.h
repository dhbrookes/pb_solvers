//
//  SHCalcUnitTest.h
//  pb_solvers_code
//
//  Created by Marielle Soniat on 10/5/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#ifndef pb_solvers_code_BesselCalcUnitTest_h
#define pb_solvers_code_BesselCalcUnitTest_h

#include "BesselCalc.h"

const int nPol = 10 ;
double precLim = 1.0e-4;



BesselConstants bConstTest = BesselConstants( nPol );
BesselCalc bCalcTest = BesselCalc( nPol, &bConstTest ) ; 


class BesselConstUTest : public ::testing::Test
{
  protected :

  double kPreFactors[10] = { -1. , 0.33333333, 0.06666667,
             0.02857143, 0.01587302, 0.01010101,
             0.00699301, 0.00512821, 0.00392157,
             0.00309598 };

  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(BesselConstUTest, first10)
{
  ASSERT_EQ( bConstTest.get_n() , nPol ); // make sure numVals stores right

  for (int i = 0; i < nPol; i++) // check that our prefactors are right
  {
    EXPECT_NEAR( bConstTest.get_kconst_val(i), kPreFactors[i], precLim);
  }
}



class BesselCalcUTest : public ::testing::Test
{
  protected :

  // for z = 1.0 and z = 10.0, calculated from python pbam_unit_test.py
  double i1[10] = {1.17520119e+00,   1.10363832e+00,   1.07344305e+00,
        1.05683451e+00,   1.04633846e+00,   1.03910889e+00,
        1.03382786e+00,   1.02980185e+00,   1.02663135e+00,
        1.02407006e+00};

  double i10[10] = {1.10132329e+03,   2.97357289e+02,   1.20594900e+02,
        6.18668362e+01,   3.69986800e+01,   2.46194746e+01,
        1.77022637e+01,   1.34885612e+01,   1.07449415e+01,
        8.86189182e+00};

  double k1[10] = {1.00000000e+00,   2.00000000e+00,   2.33333333e+00,
        2.46666667e+00,   2.53333333e+00,   2.57248677e+00,
        2.59807600e+00,   2.61606542e+00,   2.62938888e+00,
        2.63964796e+00 };

  double k10[10] = {1.00000000e+00,   1.10000000e+01,   4.43333333e+01,
        1.17666667e+02,   2.44333333e+02,   4.31105820e+02,
        6.77907167e+02,   9.79379768e+02,   1.32702447e+03,
        1.71109497e+03};
    
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(BesselCalcUTest, first10)
{
   const vector<double> mBFI10 = bCalcTest.calc_mbfI( nPol, 10.0 );
   const vector<double> mBFK10 = bCalcTest.calc_mbfK( nPol, 10.0 );
   const vector<double> mBFI1  = bCalcTest.calc_mbfI( nPol,  1.0 );
   const vector<double> mBFK1  = bCalcTest.calc_mbfK( nPol,  1.0 );

  for (int i = 0; i < nPol; i++) 
  {
    EXPECT_NEAR( mBFI10[i] , i10[i] , precLim);
    EXPECT_NEAR( mBFK10[i] , k10[i] , precLim);
    EXPECT_NEAR( mBFI1[i]  , i1[i]  , precLim);
    EXPECT_NEAR( mBFK1[i]  , k1[i]  , precLim);
  }
}

#endif

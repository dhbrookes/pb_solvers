//
//  SHCalcUnitTest.h
//  pb_solvers_code
//
//  Created by Marielle Soniat on 10/5/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#ifndef pb_solvers_code_SHCalcUnitTest_h
#define pb_solvers_code_SHCalcUnitTest_h

#include "SHCalc.h"

const int nVals = 10 ;
double preclim = 1.0e-4;

SHCalcConstants SHConstTest( nPol );


class SHCalcUTest : public ::testing::Test
{
  protected :
    
  // ensure proper calculation of (2*l-1)!!
  double doubleFactorial[10] = { 1.00000000e+00,   1.00000000e+00,
        3.00000000e+00, 1.50000000e+01,   1.05000000e+02,
        9.45000000e+02, 1.03950000e+04,   1.35135000e+05,
        2.02702500e+06, 3.44594250e+07 };
    
  // ensure propler calculation of shConsts
  double SHConst1[1]  = {};
  double SHConst5[1]  = {};
  double SHConst10[1] = {};
    
  virtual void SetUp() {}
  virtual void TearDown() {}
};


TEST_F(SHCalcUTest, first10)
{
  ASSERT_EQ( SHConstTest.get_n() , nVals ); // make sure numVals stores right

  for (int i = 0; i < nPol; i++) // check that our prefactors are right
  {
    EXPECT_NEAR( SHConstTest.get_dub_fac_val(i), doubleFactorial[i], preclim);
  }
}

#endif

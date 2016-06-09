//
//  SolverUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 6/8/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef SolverUnitTest_h
#define SolverUnitTest_h

#include "Solver.h"

/*
 Class for unit testing solver matrices
 */
class SolverUTest : public ::testing::Test
{
  protected :
  
  // ensure proper calculation of (2*l-1)!!
  double doubleFactorial[10] = { 1.00000000e+00,   1.00000000e+00,
    3.00000000e+00, 1.50000000e+01,   1.05000000e+02,
    9.45000000e+02, 1.03950000e+04,   1.35135000e+05,
    2.02702500e+06, 3.44594250e+07 };
  
  // ensure proper calculation of legendre Consts
  double LegConst1N5[6] = {1.8, 2.25,  3.0 ,  4.5 ,  9.0 , 0.0};
  double LegConst2N5[6] = {0.8, 1.25,  2.0 ,  3.5 ,  8.0 , 0.0};
  
  // ensure proper calculation of shConstant
  double SHConst[10] = {  1.00000000e+00,   1.05409255e-01,   1.12366644e-02,
    1.22602060e-03,   1.38819496e-04,   1.65921034e-05,    2.14203133e-06,
    3.09175592e-07,   5.30231766e-08,   1.24976826e-08};
  
  virtual void SetUp()
  {
    SHCalcConstants SHConstTest( nvals );
  }
  virtual void TearDown() {}
  
  SHCalcConstants SHConstTest_;
public:
  SHConstUTest( ) : SHConstTest_( nvals ) {  }
  
};


TEST_F(SolverUTest, SolvTest)
{
  ASSERT_EQ( SHConstTest_.get_n() , nvals ); // make sure numVals stores right
  
  for (int i = 0; i < nvals; i++) // check that our prefactors are right
  {
    EXPECT_NEAR( SHConstTest_.get_dub_fac_val(i), doubleFactorial[i], preclim);
  }
}

#endif /* SolverUnitTest_h */

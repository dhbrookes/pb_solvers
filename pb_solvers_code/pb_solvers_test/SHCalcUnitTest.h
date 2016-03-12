//
//  SHCalcUnitTest.h
//  pb_solvers_code
//
//  Created by Marielle Soniat on 10/5/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#ifndef SHCalcUnitTest_h
#define SHCalcUnitTest_h

#include "SHCalc.h"

/*
 Class for unit testing spherical harmonics constants
 */
class SHConstUTest : public ::testing::Test
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


TEST_F(SHConstUTest, constantTest)
{
  ASSERT_EQ( SHConstTest_.get_n() , nvals ); // make sure numVals stores right

  for (int i = 0; i < nvals; i++) // check that our prefactors are right
  {
    EXPECT_NEAR( SHConstTest_.get_dub_fac_val(i), doubleFactorial[i], preclim);
  }
}

TEST_F(SHConstUTest, LegConsTest)
{
  for (int constIt = 0; constIt < 6; constIt++)
  {
    EXPECT_NEAR( SHConstTest_.get_leg_consts1_val( 5, constIt),
                 LegConst1N5[constIt], preclim);
    EXPECT_NEAR( SHConstTest_.get_leg_consts2_val( 5, constIt),
                LegConst2N5[constIt], preclim);
  }
}


TEST_F(SHConstUTest, SHConsTest)
{
  for (int i = 0; i < nvals; i++) // check that our prefactors are right
  {
    EXPECT_NEAR( SHConstTest_.get_sh_consts_val( nvals - 1, i ),
                SHConst[i], preclim);
  }
}



/*
 Class for unit testing spherical harmonics
 */
class SHCalcUTest : public ::testing::Test
{
  protected :
  SHCalcConstants SHConstTest_;
  SHCalc SHCalcTest_;
  
public:
  SHCalcUTest( ) : SHConstTest_(10)
  {
    SHCalc SHCalcTest_(10, make_shared<SHCalcConstants> (SHConstTest_));
  }
};

TEST_F(SHCalcUTest, legendre_0)
{
  SHCalcTest_.calc_sh( 0.0, 0.0 );
  
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(        0, 0), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1, 0), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1, nvals - 1),
              0.0, preclim);
}

TEST_F(SHCalcUTest, legendre_pi3)
{
  SHCalcTest_.calc_sh( M_PI/3.0, 0.0 );
  double largeLeg = -9.44242865e+06;
  
  EXPECT_NEAR( SHCalcTest_.get_legendre_result( 0, 0), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1,        0),
              -2.67898560e-01, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1,nvals - 1)/largeLeg,
              1.0, preclim);
}

TEST_F(SHCalcUTest, legendre_2pi3)
{
  SHCalcTest_.calc_sh( 2.0*M_PI/3.0, 0.0 );
  double largeLeg = -9.44242865e+06;
  
  EXPECT_NEAR( SHCalcTest_.get_legendre_result( 0, 0), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1, 0),
              2.67898560e-01, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1,nvals - 1)/largeLeg,
              1.0, preclim);
}

TEST_F(SHCalcUTest, legendre_pi)
{
  SHCalcTest_.calc_sh( M_PI, 0.0 );
  
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(        0, 0), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(        5, 0),-1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_legendre_result(nvals - 1,nvals - 1), 0.0,
              preclim);
}

TEST_F(SHCalcUTest, sphHarm_t0p0)
{
  SHCalcTest_.calc_sh( 0.0, 0.0 );
  
  EXPECT_NEAR( SHCalcTest_.get_result( 0, 0).real(), 1.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 0, 0).imag(), 0.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 5, 5).real(), 0.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 5, 5).imag(), 0.0, preclim);
}


TEST_F(SHCalcUTest, sphHarm_t05p05)
{
  SHCalcTest_.calc_sh( 0.5, 0.5 );
  
  EXPECT_NEAR( SHCalcTest_.get_result( 5, 0).real(),-0.16928726, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 0, 0).imag(), 0.0, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 5, 5).real(),-0.010066221, preclim);
  EXPECT_NEAR( SHCalcTest_.get_result( 5, 5).imag(), 0.00751969, preclim);
}



#endif

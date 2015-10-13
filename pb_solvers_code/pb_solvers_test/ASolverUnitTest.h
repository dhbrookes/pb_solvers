//
//  ASolverUnitTest.h
//  pbsolvers
//
//  Created by Marielle Soniat on 10/6/15.
//  Copyright (c) 2015 Marielle Soniat. All rights reserved.
//

#ifndef pbsolvers_ASolverUnitTest_h
#define pbsolvers_ASolverUnitTest_h

#include <iostream>
#include "ASolver.h"
using namespace std;

class ASolverUTest : public ::testing::Test
{
public :
  
protected :
  
  int vals_;
  Constants const_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    EPt pos[2]   = { EPt( 0.0, 0.0, -5.0 ), EPt( 10.0, 7.8, 25.0 ) };
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(1);
      vector<EPt> posCharges(1);
      
      charges[0] = cg[molInd];
      posCharges[0] = pos[molInd];
      
      Molecule molNew = Molecule( M, rd[molInd], charges, posCharges );
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
} ;

// functions to test:

TEST_F(ASolverUTest, checkGamma)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( vals );
  BesselCalc bCalcu        = BesselCalc( vals, &bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( vals );
  SHCalc SHCalcu           = SHCalc( vals, &SHConsta );
  System sys               = System( const_, mol_ );
  
  ASolver ASolvTest        = ASolver( nmol, vals, &bCalcu, &SHCalcu, sys );

  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 0),  1.463995711, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 4),  1.760111936, preclim);
  
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 1),  1.621243794, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 6),  1.799701878, preclim);
}

TEST_F(ASolverUTest, checkDelta)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( vals );
  BesselCalc bCalcu        = BesselCalc( vals, &bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( vals );
  SHCalc SHCalcu           = SHCalc( vals, &SHConsta );
  System sys               = System( const_, mol_ );
  
  ASolver ASolvTest        = ASolver( nmol, vals, &bCalcu, &SHCalcu, sys );
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 0)/56.03476045, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 4)/73361234.99, 1.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 1)/46846.22401, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 6)/8.00377E+14, 1.0, preclim);
}
 

TEST_F(ASolverUTest, checkE)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( vals );
  BesselCalc bCalcu        = BesselCalc( vals, &bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( vals );
  SHCalc SHCalcu           = SHCalc( vals, &SHConsta );
  System sys               = System( const_, mol_ );
  
  ASolver ASolvTest        = ASolver( nmol, vals, &bCalcu, &SHCalcu, sys );
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 0, 0).real(), 5.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 0, 0).imag(), 0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 5,5+0).real()/-15625, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 5,5+0).imag(),        0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 6,6-5).real(),        0.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 6,6-5).imag(),        0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 0, 0).real(),-0.4, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 0, 0).imag(), 0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 3,3-3).real()/184.52,  1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 3,3-3).imag()/417.127, 1.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 6,6-5).real()/5.31968e+06, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 6,6-5).imag()/-916110,     1.0, preclim);

}

/*
 
 // testing matrix assertions for accessing matrix members
 TEST_F(MyMatrixUTest, exceptAccess)
 {
 MyMatrix <double> testMat_( testInp_ );

 ASSERT_THROW(testMat_( 0,-1), MatrixAccessException);
 }
 
 
 */



#endif

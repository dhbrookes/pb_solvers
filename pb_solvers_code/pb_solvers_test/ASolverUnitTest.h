//
//  ASolverUnitTest.h
//  pbsolvers
//
//  Created by Marielle Soniat on 10/6/15.
//  Copyright (c) 2015 Marielle Soniat. All rights reserved.
//

#ifndef pbsolvers_ASolverUnitTest_h
#define pbsolvers_ASolverUnitTest_h

#include "ASolver.h"

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
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.5 ), Pt( 11.0, 6.9, 24.3 ) };
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(1);
      vector<Pt> posCharges(1);
      
      charges[0]    = cg[molInd];
      posCharges[0] = cgPos[molInd];
      
      Molecule molNew( M, rd[molInd], charges, posCharges, pos[molInd]);
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
} ; // end ASolverUTest


TEST_F(ASolverUTest, checkGamma)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);

  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 1).real(),  1.463995711, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 5).real(),  1.760111936, preclim);
  
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 2).real(),  1.621243794, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 7).real(),  1.799701878, preclim);
}

TEST_F(ASolverUTest, checkDelta)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 1).real()/56.03476045, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 5).real()/73361234.99, 1.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 2).real()/46846.22401, 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 7).real()/8.00377E+14, 1.0, preclim);
}
 

TEST_F(ASolverUTest, checkE)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 0, 0).real(), 5.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 5, 0).real(),   -0.15625, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 5, 0).imag(),        0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 6, -5).real(),        0.0, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 0, 6, -5).imag(),        0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 0, 0).real(),-0.4, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 3, -2).real(), 0.0728481807168, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 3, -2).imag(), 0.6901393135990, preclim);

  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 6, -5).real(), -1.75961914197, preclim);
  EXPECT_NEAR( ASolvTest.get_E_ni( 1, 6, -5).imag(), -1.01329895122, preclim);
}

TEST_F(ASolverUTest, checkSH)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 0, 0).real(), 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 5, 0).real(),-1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 5, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 6,-5).real(), 0.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 6,-5).imag(), 0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 0, 0).real(), 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 3,-2).real(),-0.0522110883, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 3,-2).imag(), 0.4946303982, preclim);
  
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 6, 5).real(), 0.3615486465, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 6, 5).imag(),-0.2082023636, preclim);
  
}

TEST_F(ASolverUTest, checkT)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  
}

TEST_F(ASolverUTest, checkA)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(20);
}

#endif

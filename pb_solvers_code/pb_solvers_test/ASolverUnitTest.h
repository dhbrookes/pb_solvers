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
  vector< Molecule > mol_sing_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.5 ), Pt( 11.0, 6.9, 24.3 ) };
    
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, -0.0 ) };
    
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
      
      charges[0]    = 2.0;
      posCharges[0] = cgPosSi[molInd];

      Molecule molSing( M, 10.0, charges, posCharges);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  double A0[15] = {5.06492124, -3.6629753, -0.000641312455, 2.01744339,
    -0.000150832137,-5.86932378e-06,-1.05487251,-1.46904154e-05,8.44898379e-07,
    8.34955966e-07, 0.541011058, 4.98503746e-07, 9.23014214e-07, 1.77640914e-07,
    3.51977104e-09,};
  
  double A0_im[15] = { -4.51751614e-09, -0.000253764627, -0.000546361723,
    -0.000140256989, -0.000161954107, -4.28036153e-05, -3.85518402e-05,
    -2.87295411e-05, -8.62026893e-06, -9.90550946e-07, -7.87063036e-06,
    -4.16328095e-06, -5.59013e-07, 8.23515067e-08, 2.35969719e-08};
  
   double A1[15] = { -0.402210927, 0.272084187, -0.448632042, 0.332513859,
     0.584393585, -0.0735997277, -1.07437691, -0.058387728, 0.121844253,
     0.542085518, 0.832263495, -0.937549973, -0.0580719897, -1.02883583,
     1.16366737};
  
  double A1_im[15] = {-6.00898647e-07,1.0270678e-05,0.348521971,-9.72430906e-06,
    -0.478199881, 0.722212407, 6.23391427e-06, 0.0283916673, -1.17520566,
    0.745545848, -3.18354738e-06, 0.853352497, 0.562026511, -1.41661515,
    0.248528357};

//  double A0[15] = { 5.06496536, -3.6619316, -0.000503250061, 2.01779987,
//    -9.82900676e-05, -1.0682227e-05, -1.05479692, -3.28880159e-06,
//    -2.71180395e-06, -4.4780814e-07, 0.541023412, 2.40308854e-06,
//    -3.44083312e-07, -4.12058842e-07, -8.58237216e-08};
//  
//  double A0_im[15] = {0, 0.000340562962, -0.000330632776, 0.000163525329,
//    -2.30829341e-05, -1.86457723e-05, 3.5928091e-05, 1.53052678e-05,
//    2.2200138e-06, 9.73474763e-09, 5.13843793e-06, 5.61706783e-06,
//    2.48309863e-06, 4.65173384e-07, 2.04983496e-08 };
//  
//  double A1[15] = {-0.401833604, 0.265698439, -0.443888271, 0.338291377,
//    0.58064118,-0.07417005,-1.0778658, -0.0567916512, 0.122326234,
//    0.54187217, 0.833960226,-0.9379799,-0.0583482748,-1.02863481, 1.16375327};
//  
//  double A1_im[15] = {0, 0, 0.352227876, 0, -0.481130826, 0.71993116, 0,
//    0.0296375441, -1.1732774, 0.746026633, 0, 0.85301775, 0.560921204,
//    -1.417068, 0.24848249 };
  
  double A0Sing[15] = {2.08493611, 0.0595090739,0, 0.026460114, 0, 0,
    0.00906667786, 0, 0, 0, 0.00287090794, 0, 0, 0, 0};
  
  double A1Sing[15] = {2.08493623, -0.059510451, 0, 0.026460955, 0, 0,
    -0.00906706693, 0, 0, 0, 0.00287106903, 0, 0, 0, 0};

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
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 1).real(),   0.87554313, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 5).real(),   0.06832297, preclim);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 2).real(),   11.4370663, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 7).real()/181.9847, 1.0, preclim);
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
  ASolvTest.solve_A(1E-20);
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).real(), A0[ct],    preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).imag(), A0_im[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).real(), A1[ct],    preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).imag(), A1_im[ct], preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkASing)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-7);
  
  int ct = 0;

  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).real(), A0Sing[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).imag(),        0.0, preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).real(), A1Sing[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).imag(),        0.0, preclim);
      ct++;
    }
  }
}

#endif

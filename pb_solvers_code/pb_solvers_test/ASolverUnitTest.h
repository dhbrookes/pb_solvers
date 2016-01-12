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
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, 0.0 ) };
    
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
  
//  double A0[15] = {5.06492124, -3.6629753, -0.000641312455, 2.01744339,
//    -0.000150832137,-5.86932378e-06,-1.05487251,-1.46904154e-05,8.44898379e-07,
//    8.34955966e-07, 0.541011058, 4.98503746e-07, 9.23014214e-07, 1.77640914e-07,
//    3.51977104e-09,};
//  double A0_im[15] = { -4.51751614e-09, -0.000253764627, -0.000546361723,
//    -0.000140256989, -0.000161954107, -4.28036153e-05, -3.85518402e-05,
//    -2.87295411e-05, -8.62026893e-06, -9.90550946e-07, -7.87063036e-06,
//    -4.16328095e-06, -5.59013e-07, 8.23515067e-08, 2.35969719e-08};
//   double A1[15] = { -0.402210927, 0.272084187, -0.448632042, 0.332513859,
//     0.584393585, -0.0735997277, -1.07437691, -0.058387728, 0.121844253,
//     0.542085518, 0.832263495, -0.937549973, -0.0580719897, -1.02883583,
//     1.16366737};
//  double A1_im[15] = {-6.00898647e-07,1.0270678e-05,0.348521971,-9.72430906e-06,
//    -0.478199881, 0.722212407, 6.23391427e-06, 0.0283916673, -1.17520566,
//    0.745545848, -3.18354738e-06, 0.853352497, 0.562026511, -1.41661515,
//    0.248528357};

  double A0[15] = { 5.06496536, -3.6619316, -0.000503250061, 2.01779987,
    -9.82900676e-05, -1.0682227e-05, -1.05479692, -3.28880159e-06,
    -2.71180395e-06, -4.4780814e-07, 0.541023412, 2.40308854e-06,
    -3.44083312e-07, -4.12058842e-07, -8.58237216e-08};
  double A0_im[15] = {0, 0.000340562962, -0.000330632776, 0.000163525329,
    -2.30829341e-05, -1.86457723e-05, 3.5928091e-05, 1.53052678e-05,
    2.2200138e-06, 9.73474763e-09, 5.13843793e-06, 5.61706783e-06,
    2.48309863e-06, 4.65173384e-07, 2.04983496e-08 };
  double A1[15] = {-0.401833604, 0.265700049, -0.449876813, 0.338288496,
    0.586157705,-0.0734981144, -1.07786325, -0.0598377042, 0.12171614,
    0.54211679, 0.83395858, -0.93662999, -0.057964058, -1.028876, 1.1636555};
  double A1_im[15] = {0, 4.58511161e-06, 0.347557646, -8.25770523e-06,
    -0.47683061, 0.722606727, 7.39004703e-06, 0.0272651594, -1.17570528,
    0.74547344,-4.85359351e-06, 0.85406746, 0.56244862,-1.4165226, 0.24853490};
  
  double A0Sing[15] = {2.08493611, 0.0595090739,0, 0.026460114, 0, 0,
    0.00906667786, 0, 0, 0, 0.00287090794, 0, 0, 0, 0};
  double A1Sing[15] = {2.08493623, -0.059510451, 0, 0.026460955, 0, 0,
    -0.00906706693, 0, 0, 0, 0.00287106903, 0, 0, 0, 0};
  
  double dA00[15] = {1.0544637e-06,7.85368589e-05,-3.7739641e-05,1.1710611e-05,
    -4.2581267e-06,-2.54957295e-06,1.1509835e-06,-2.02804679e-07,-3.6652067e-07,
    -4.9837825e-08,8.9749664e-08,2.9501223e-09,-3.6136045e-08,-8.8566072e-09,0};
  double dA01im[15] = {0,0,4.49886716e-05,0,5.90910631e-06,1.96972506e-06,0,
    4.233064e-07,2.48485723e-07,5.45274128e-09,0,1.98818687e-08,1.98691187e-08,
    0,-2.28180158e-09};
  double dA11[15] = {-3.1633912e-06,0.00015605884,5.5533741e-05,-1.438732e-05,
    -9.82898223e-06,-8.25866317e-07,8.00828799e-07,1.0218887e-06,
    1.40543141e-07,-4.82103708e-08,-2.17563271e-08,-8.19726265e-08,
    -1.62881095e-08,8.72542217e-09,3.2623144e-09};
  double dA10im[15] = {0,0,-1.44387654e-05,0,3.28980681e-06,-1.10667758e-06,0,
    -4.39392408e-07,1.06332252e-07,9.822109e-08,0,4.54976027e-08,-3.4394996e-09,
    -1.4438813e-08,-3.54879966e-09};
  
  double dATrip02[15] = {2.5859665e-05,0.000620232366,0.000285436615,
    0.000236529621,8.34048199e-05,-4.91547478e-05,6.02966515e-05,1.00209838e-05,
    -1.8599944e-05,2.4739817e-06,1.2321061e-05,-5.74248967e-07,-4.85496063e-06,
    1.14977575e-06,9.95493175e-08};
  double dATrip11im[15] = {1.26068667e-07,1.77759087e-06,-0.0183902725,
    3.50966837e-06,0.00404630668,0.031857372,-4.92548117e-06,0.00782664265,
    -0.00183388625,-0.0429967544,5.87124061e-06,0.000568727219,-0.018891895,
    -7.10651767e-05,0.051565011};
  double dATrip21[15] = {1.11111691e-06,2.36981299e-05,-2.85249661e-05,
    -2.31735024e-06,2.95909727e-06,2.85184934e-06,-1.8419939e-07,4.31283004e-08,
    1.5601647e-07,3.89577124e-08,4.5967117e-09,-1.62540654e-08,-3.85136191e-09,
    1.48585811e-08,4.41682557e-09};
  double dATrip10im[15] = {-3.62501039e-08,-3.17816832e-07,-0.000808271535,
    -1.97785916e-07,0.00276094976,-0.000692074481,3.07878903e-06,-0.00162795109,
    0.00110762117,0.000646125447,-2.84638226e-06,0.00264295314,0.000193963401,
    -0.000577338053,-0.000156288848};
  
  double dASing00[15] = {0,0,-0.0012172202,0,-0.0009384587,0,0,-0.00045537497,
    0,0,0,-0.000186430392,0,0,0};
  double dASing00im[15] = {0,0,-0.00121722018,0,-0.00093845868,0,0,
    -0.000455374968,0,0,0,-0.000186430392,0,0,0};
  double dASing01[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double dASing01im[15] = {0,0,0.00121722018,0,0.00093845868,0,0,
    0.000455374968,0,0,0,0.000186430392,0,0,0};
  double dASing11[15] = {-0.000411266015,0.0044280306,0,-0.00261408043,0,0,
    0.00113979399,0,0,0,-0.000444103676,0,0,0,0};
  double dASing11im[15] = {0,0,-0.00121722018,0,0.00093845868,0,0,
    -0.000455374968,0,0,0,0.000186430392,0,0,0};

  double L0[15] = {-0.0042958084,-0.0015199557,-0.00039314106,-0.000151941346,
    -0.0001196473,-1.2948965e-05,3.9843662e-05,-7.88050185e-06,-6.0250525e-06,
    -9.78873154e-07,2.38774686e-05,9.83466712e-06,-1.48977456e-06,
    -1.73451619e-06,-3.6075022e-07};
  double L0_im[15] = {0,0.00026570391,-0.00025835448,0.000197697415,
    -2.85442258e-05,-2.2679504e-05,7.92585234e-05,3.32760866e-05,4.72712042e-06,
    1.12487091e-09,2.17628492e-05,2.35222037e-05,1.03871094e-05,1.94603771e-06,
    8.53505184e-08};
  double L1[15] = {0.0413679257,-0.0173166068,-0.00405117907,0.003729965,
    0.00162800449,0.000105438302,-0.00064643838,-0.00045355362,-4.81366529e-05,
    1.3683007e-05,9.4098892e-05,0.000105662438,1.58220788e-05,-7.0246614e-06,
    -2.26791712e-06};
  double L1_im[15] = {0,5.44303496e-07,-0.00315982074,-4.45344828e-07,
    0.00126969936,0.000419966072,2.1175274e-07,-0.000353670484,-0.000191687314,
    -3.09373544e-05,-7.75938102e-08,8.23654847e-05,6.2979024e-05,1.58859638e-05,
    1.2155081e-06};
  
  double L0Sing[15] = {0.0210272529,0.0124511487,0,0.00424544438,0,0,
    0.0013056269,0,0,0,0.000389723882,0,0,0,0};
  double L0SingIm[15] = {};
  double L1Sing[15] = {0.0210272529,-0.0124511487,0,0.00424544438,0,0,
    -0.0013056269,0,0,0,0.000389723882,0,0,0,0};
  double L1SingIm[15] = {};
  
  double dLdx0[15] = {1.6750432e-06,-6.48184643e-05,2.34002566e-05,
    9.24707607e-06,1.70674585e-05,4.63747276e-06,2.14310222e-06,4.22726785e-06,
    1.31460586e-06,1.10177468e-07,-4.90918167e-07,-2.74985149e-06,
    -8.99415319e-07,1.44385542e-07,1.00536605e-07};
  double dAdy0im[15] = {-1.61770896e-06,-3.41962773e-06,3.5234244e-05,
    -1.05669019e-05,7.2297967e-06,4.09008335e-06,-3.47883114e-06,7.54460728e-07,
    1.96413284e-06,4.91172588e-07,-1.26891901e-08,-2.75940448e-06,
    -1.43810291e-06,-9.21941465e-08,6.53481772e-08};
  double dLdz0[15] = {-0.000150691983,-1.33519347e-05,-2.06166021e-05,
    5.03159396e-05,2.06163519e-05,1.07343682e-06,2.0450875e-05,1.79780738e-05,
    2.12705147e-06,-5.04752566e-07,-4.9965286e-06,-4.07149703e-06,
    -6.09290158e-07,1.72960474e-07,5.43921381e-08};
  double dLdx1[15] = {-0.000987674279,0.000187867261,-0.000357770967,
    -0.000293354883,5.61544215e-05,4.6981816e-05,3.09724418e-05,-3.49290618e-05,
    -2.36043496e-05,-2.24351551e-06,-2.44654805e-05,-6.57150526e-06,
    5.77399697e-06,1.98995506e-06,-1.03763066e-08};
  double dAdy1im[15] = {2.61723492e-07,-3.77632927e-08,-0.00037629043,
    9.4750488e-08,9.7902191e-05,3.25167845e-05,-1.47164999e-08,-3.90081e-05,
    -2.15522995e-05,-1.42788374e-06,1.30984974e-08,1.38167558e-07,
    9.21538861e-07,-7.45070972e-07,-4.6128438e-07};
  double dAdz1im[15] = {-1.18032719e-07,-7.74616576e-08,0.000240731667,
    7.46244492e-08,-0.000116771598,-5.26094072e-05,-4.6568393e-08,
    4.78806389e-05,3.23462883e-05,6.20808799e-06,1.75724926e-08,-7.2420955e-06,
    -8.11121113e-06,-2.60581756e-06,-2.39864327e-07};
  
  double dLdx0Sing[15] = {0,0,-0.00025467424,0,-0.000150567895,0,0,
    -6.5572367e-05,0,0,0,-2.53063413e-05,0,0,0};
  double dAdy0imSing[15] = {0,0,-0.00025467424,0,-0.000150567895,0,0,
    -6.5572367e-05,0,0,0,-2.53063413e-05,0,0,0};
  double dLdz0Sing[15] = {0.00128088088,0.000940448055,0,0.000393219492,0,0,
    0.000151918661,0,0,0,5.43212252e-05,0,0,0,0};
  double dLdx1Sing[15] = {0,0,-0.000254674238,0,0.000150567895,0,0,
    -6.5572367e-05,0,0,0,2.53063413e-05,0,0,0};
  double dAdy1imSing[15] = {0,0,-0.00025467424,0,0.00015056789,0,0,
    -6.5572367e-05,0,0,0,2.53063413e-05,0,0,0};
  double dAdz1imSing[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  
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
  ASolvTest.solve_A(1E-20);
  
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

TEST_F(ASolverUTest, checkgradA)
{
  vector<double> charges(1); vector<Pt> posCharges(1);
  charges[0] = 2.0; posCharges[0] = Pt(-10.0, 7.8, 25.0);
  Molecule molNew( 1, 2.0, charges, posCharges, posCharges[0]);
  mol_.push_back( molNew );
  
  const int vals           = nvals;
  int nmol                 = 3;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-16);

  int ct = 0;
  for ( int n = 0; n < 3; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      if (dATrip02[ct] != 0)
        EXPECT_NEAR( ASolvTest.get_dAdx_ni( 0, 2, n, m).real()/dATrip02[ct],
                     1.0, preclim);
      if (dATrip11im[ct] != 0)
        EXPECT_NEAR( ASolvTest.get_dAdy_ni( 1, 1, n, m).imag()/dATrip11im[ct],
                     1.0, preclim);
      if (dATrip21[ct] != 0)
        EXPECT_NEAR( ASolvTest.get_dAdz_ni( 2, 1, n, m).real()/dATrip21[ct],
                     1.0, preclim);
      if (dATrip10im[ct] != 0)
        EXPECT_NEAR( ASolvTest.get_dAdx_ni( 1, 0, n, m).imag()/dATrip10im[ct],
                     1.0, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkgradASing)
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
  ASolvTest.solve_A(1E-20);
  ASolvTest.solve_gradA(1E-16);
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( ASolvTest.get_dAdx_ni( 0, 0, n, m).real(),
                  dASing00[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 0, n, m).imag(),
                  dASing00im[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 1, n, m).real(),
                  dASing01[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 1, n, m).imag(),
                  dASing01im[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_dAdz_ni( 1, 1, n, m).real(),
                  dASing11[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 1, 1, n, m).imag(),
                  dASing11im[ct], preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkL)
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
  VecOfMats<cmplx>::type myL = ASolvTest.calc_L();
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      if (L0[ct] != 0)
        EXPECT_NEAR( myL[0](n,m+nvals).real()/L0[ct],    1.0, preclim);
      if (L0_im[ct] != 0)
        EXPECT_NEAR( myL[0](n,m+nvals).imag()/L0_im[ct], 1.0, preclim);
      if (L1[ct] != 0)
        EXPECT_NEAR( myL[1](n,m+nvals).real()/L1[ct],    1.0, preclim);
      if (L1_im[ct] != 0)
        EXPECT_NEAR( myL[1](n,m+nvals).imag()/L1_im[ct], 1.0, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkLSing)
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
  ASolvTest.solve_A(1E-20);
  
  VecOfMats<cmplx>::type myL = ASolvTest.calc_L();
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( myL[0](n,m+nvals).real(), L0Sing[ct], preclim);
      EXPECT_NEAR( myL[0](n,m+nvals).imag(),          0, preclim);
      EXPECT_NEAR( myL[1](n,m+nvals).real(), L1Sing[ct], preclim);
      EXPECT_NEAR( myL[1](n,m+nvals).imag(),          0, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkdL)
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
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  MyVector<VecOfMats<cmplx>::type > mydL = ASolvTest.calc_gradL();
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      if (dLdx0[ct] != 0)
        EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0[ct],   preclim);
      if (dAdy0im[ct] != 0)
        EXPECT_NEAR( mydL[0][1](n,m+nvals).imag(), dAdy0im[ct], preclim);
      if (dLdx0[ct] != 0)
        EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0[ct],   preclim);
      if (dLdx1[ct] != 0)
        EXPECT_NEAR( mydL[1][0](n,m+nvals).real(), dLdx1[ct],   preclim);
      if (dAdy1im[ct] != 0)
        EXPECT_NEAR( mydL[1][1](n,m+nvals).imag(), dAdy1im[ct], preclim);
      if (dAdz1im[ct] != 0)
        EXPECT_NEAR( mydL[1][2](n,m+nvals).imag(), dAdz1im[ct], preclim);
      ct++;
    }
  }  
}

TEST_F(ASolverUTest, checkdLSing)
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
  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
  
  MyVector<VecOfMats<cmplx>::type > mydL = ASolvTest.calc_gradL();
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0Sing[ct],   preclim);
      EXPECT_NEAR( mydL[0][1](n,m+nvals).imag(), dAdy0imSing[ct], preclim);
      EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0Sing[ct],   preclim);
      EXPECT_NEAR( mydL[1][0](n,m+nvals).real(), dLdx1Sing[ct],   preclim);
      EXPECT_NEAR( mydL[1][1](n,m+nvals).imag(), dAdy1imSing[ct], preclim);
      EXPECT_NEAR( mydL[1][2](n,m+nvals).imag(), dAdz1imSing[ct], preclim);
      ct++;
    }
  }
}

#endif

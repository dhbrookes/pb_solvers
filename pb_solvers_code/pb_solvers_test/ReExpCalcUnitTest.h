//
//  ReExpCalcUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/13/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ReExpCalcUnitTest_h
#define ReExpCalcUnitTest_h

#include <iostream>
#include "ReExpCalc.h"

class ReExpUTest : public ::testing::Test
{
public :
  
  ReExpUTest() : vals_(nvals)   {  }
  
  int get_vals()                { return vals_; }
  
protected :
  
  int vals_;
  
  double R0Zreal[nvals]  = { 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double R0Zimag         =   0.0 ;
  
  double RM5Zreal[nvals] = { 0, 0, 0, 0, 0,1.0, 0, 0, 0, 0};
  double RM5Zimag        =   0.0 ;
  
  double R0real[nvals] = { -0.05998438815,  0.21425438524,  0.02614110987,
    0.02614098352,  0.03502736468,  -0.2576350453,  -0.0520278953,
    0.22443080465,  0.01126460251,  0.12826498127};
  double R0imag[nvals] = { 0.0,  0.13352084876,  0.05326969298,
  -0.257709011789, -0.045282977706, 0.09554032805, -0.010664740507,
   0.21308689852,   0.04338655884, -0.409892178514 };
  
  double RM1Real[nvals] = {-0.252399,  -0.027798,  0.115993,  0.002635,
    0.165055,  -0.090395,  -0.263908,  0.139210,  0.064392,  0.116713 };
  double RM7Imag[nvals] = {-0.000000,  -0.108268,  0.037856,  0.212323,
    -0.367881,  0.186130,  0.056455,  0.188322,  -0.910017,  1.550056 };
  
  double SN0Z[nvals] = {0.970147, 4.99775, 24.99617, 124.988519, 624.958996,
    3124.840538,  15624.347650,  78122.240044,  390613.040154, 1953072.235847};
  
  double SN0[nvals] = { 0.096100,  0.368228,  1.150332,  3.550246,  10.932428,
    33.638947,  103.470482,  318.206275,  978.479829,  3008.591425 };
  
  double SNNZ[nvals] = {0.197923,  0.000000,  -0.000764,  0.001274,  -0.002677,
    0.006426,  -0.016831,  0.046891,  -0.136772,  0.413374  };
  
  double SNN[nvals] = {0.368228, 0,-0.339689,5.36844,-106.599514,2418.548197,
    -59872.902182,  1576398.006538,  -43454097.266868,  1241141220.045197 };
  
  double SMMZ[nvals] = {0.19923861,  0.48877513,  0.77332760,  1.05695152,
    1.34024129, 1.62337339, 1.90641867, 2.18941109, 2.47236893, 0.00000000 };

  double SMM[nvals] = {1.15033201,  26.73265165,  400.10876211,  5170.65006781,
    61980.40129512,  709610.89326185,  7876260.175968,  85488801.789434,
    912350693.98094785, 0.0 };
  
  double SMLZ[nvals] = { 5.41735113,  10.98949976,  19.03910544,  29.85272945,
    43.68976929,  60.78900127,  0.00119130,  0.0, 0.0, 0.0 };

  double SML[nvals] = { 8626.68751768,  165395.13830094,  2708112.63388879,
    40130100.04533036,  555039663.08194208,  7298284604.20387840,  3575253.29395807,  0, 0, 0 };
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};


TEST_F(ReExpUTest, checkR0Zpt)
{
  Constants Cst;
  ShPt testPt = EPt( 0.0, 0.0, 5.0).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCo,
                           kap, lambda );

  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).real(),
                R0Zreal[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).imag(),
                R0Zimag, preclim);
  }
  
}

TEST_F(ReExpUTest, checkR0)
{
  Constants Cst;
  ShPt testPt = EPt( 6.9,-4.3,-0.2).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).real(),
                R0real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).imag(),
                R0imag[s], preclim);
  }
  
}

TEST_F(ReExpUTest, checkRZpt)
{
  Constants Cst;
  ShPt testPt = EPt( 0.0, 0.0, 5.0).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 5, s).real(),
                RM5Zreal[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 1, s).imag(),
                RM5Zimag, preclim);
  }
  
}

TEST_F(ReExpUTest, checkR)
{
  Constants Cst;
  ShPt testPt = EPt( 6.9,-4.3,-0.2).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 1, s).real(),
                RM1Real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 7, s).imag(),
                RM7Imag[s], preclim);
  }
  
}


TEST_F(ReExpUTest, checkS0Zpt)
{
  Constants Cst;
  ShPt testPt = EPt( 0.0, 0.0, 1.0).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCo,
                           kap, lambda );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_sval( 0, s, 0), SN0Z[s],               preclim);
    EXPECT_NEAR( ReExpTest.get_sval( s, 0, 0), SN0Z[s] * pow(-1.0,s), preclim);
  }
}


TEST_F(ReExpUTest, checkS0)
{
  Constants Cst;
  ShPt testPt = EPt( 6.9,-4.3,-0.2).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_sval( 0, s, 0), SN0[s],               preclim);
    EXPECT_NEAR( ReExpTest.get_sval( s, 0, 0), SN0[s] * pow(-1.0,s), preclim);
  }
  
}


TEST_F(ReExpUTest, checkSZpt)
{
  Constants Cst;
  ShPt testPt = EPt( 0.0, 0.0, 5.0).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_sval(   s,  s+1,   0), SNNZ[s], preclim);
    EXPECT_NEAR( ReExpTest.get_sval(   s,  s+2,   s), SMMZ[s], preclim);
    EXPECT_NEAR( ReExpTest.get_sval( s+3,  s+4, s+1), SMLZ[s], preclim);
  }
  
}


 
TEST_F(ReExpUTest, checkS)
{
  Constants Cst;
  ShPt testPt = EPt( 6.9,-4.3,-0.2).convert_to_spherical();
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, &ReExpCoeff,
                           kap, lambda );
  
  for ( int s = 0; s < nvals; s++ )
  {
    if ( SNN[s] != 0)
      EXPECT_NEAR( ReExpTest.get_sval( s, s+1, 0)/SNN[s], 1.0, preclim);
    
    if ( SMM[s] != 0)
      EXPECT_NEAR( ReExpTest.get_sval( s, s+2, s)/SMM[s], 1.0, preclim);
    
    if ( SML[s] != 0)
      EXPECT_NEAR( ReExpTest.get_sval( s+3, s+4, s+1)/SML[s], 1.0, preclim);
  }
  
}



#endif /* ReExpCalcUnitTest_h */

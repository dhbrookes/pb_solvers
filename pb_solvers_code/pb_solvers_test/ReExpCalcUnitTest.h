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
  
  double SN0Zreal[nvals] = {0.970147, 4.99775, 24.99617, 124.988519, 624.958996,
    3124.840538,  15624.347650,  78122.240044,  390613.040154, 1953072.235847};
  
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
    EXPECT_NEAR( ReExpTest.get_sval( 0, s, 0),SN0Zreal[s],             preclim);
    EXPECT_NEAR( ReExpTest.get_sval( s, 0, 0),SN0Zreal[s]*pow(-1.0,s), preclim);
  }
}

/*
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
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 5, s).real(),
                RM5Zreal[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 1, s).imag(),
                RM5Zimag, preclim);
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

*/



#endif /* ReExpCalcUnitTest_h */

//
//  ReExpCalcUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/13/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ReExpCalcUnitTest_h
#define ReExpCalcUnitTest_h

#include "ReExpCalc.h"

class ReExpUTest : public ::testing::Test
{
public :
  
  ReExpUTest() : vals_(nvals)   {  }
  
  int get_vals()                { return vals_; }
  
protected :
  
  int vals_;
  
  double R0Zreal[nvals] = { 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double R0Zimag        =   0.0 ;
  
  double R0real[nvals] = { -0.05998438815,  0.21425438524,  0.02614110987,
    0.02614098352,  0.03502736468,  -0.2576350453,  -0.0520278953,
    0.22443080465,  0.01126460251,  0.12826498127};
  double R0imag[nvals] = { 0.0,  0.13352084876,  0.05326969298,
  -0.257709011789, -0.045282977706, 0.09554032805, -0.010664740507,
   0.21308689852,   0.04338655884, -0.409892178514 };
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};


TEST_F(ReExpUTest, checkR0Zpt)
{
  Constants Cst;
  ShPt testPt = EPt( 0.0, 0.0, 5.0).convert_to_spherical();
  SHCalcConstants shCon( nvals );
  SHCalc shCalc( nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( nvals );
  BesselCalc      bCal( nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, kap, lambda );

  for ( int s = -nvals+1; s <= nvals-1; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).real(),
                R0Zreal[s+nvals-1], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).imag(),
                R0Zimag, preclim);
  }
  
}

TEST_F(ReExpUTest, checkR0)
{
  Constants Cst;
  ShPt testPt = EPt( 6.9,-4.3,-0.2).convert_to_spherical();
  SHCalcConstants shCon( nvals );
  SHCalc shCalc( nvals, &shCon );
  shCalc.calc_sh( testPt.get_theta(), testPt.get_phi());
  
  BesselConstants bCon( nvals );
  BesselCalc      bCal( nvals, &bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  
  ReExpCoeffs_IJ ReExpTest( nvals, testPt, &shMat, &bCal, kap, lambda );
  
  for ( int s = -nvals+1; s <= nvals-1; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).real(),
                R0real[s+nvals-1], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 0, s).imag(),
                R0imag[s+nvals-1], preclim);
  }
  
}


#endif /* ReExpCalcUnitTest_h */

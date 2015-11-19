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

/*
 Class to test R and S constants : A, B, alpha, beta, nu and mu
 */
class ReExpConstUTest : public ::testing::Test
{
  public :
  
  ReExpConstUTest() : vals_(nvals), kappa_(0.5), lambda_(1.0),
  ReExpCo_( kappa_, lambda_, nvals){}
  
  int get_vals()                { return vals_; }
  void set_vals( double kappa, double lambd )
  {
    kappa_  = kappa;
    lambda_ = lambd;
    ReExpCo_ = ReExpCoeffsConstants( kappa_, lambda_, vals_);
  }
  
  protected :
  int vals_;
  double kappa_;
  double lambda_;
  ReExpCoeffsConstants ReExpCo_;
  
  double A9[nvals*2] = { 0, 0.21821789,  0.30037570,  0.35751860,  0.40050094,
    0.43355498, 0.45883147, 0.47756693,  0.49051147,  0.49811675,  0.50062617,
    0.49811675, 0.49051147, 0.47756693,  0.45883147,  0.43355498,  0.40050094,
    0.35751860,  0.30037570,  0.21821789};
  
  double B9[nvals*2] = { 0,-0.97332853,-0.91766294, -0.86199423,-0.80632177,
    -0.75064472,-0.69496197,-0.63927203,-0.58357285,-0.52786151, 0.47213369,
    0.41638277, 0.36059806, 0.30476098, 0.24883630, 0.19274777, 0.13629326,
    0.07868895, 0.00000000, 0.00000000};
  
  double AL9[2*nvals] = { 0,4.35889894, 6., 7.14142843, 8., 8.66025404,
    9.16515139, 9.53939201, 9.79795897, 9.94987437, 10., 9.94987437,
    9.79795897, 9.53939201, 9.16515139, 8.66025404,  8., 7.14142843,
    6.,   4.35889894};
  
  double BEOrg9[2*nvals] = { 0., 0.00273114, 0.00375940, 0.00447458, 0.00501253,
    0.00542622, 0.00574258, 0.00597706, 0.00613907, 0.00623426, 0.00626566,
    0.00623426, 0.00613907, 0.00597706, 0.00574258, 0.00542622, 0.00501253,
    0.00447458, 0.00375940, 0.00273114 };
  
  double BE9[2*nvals] = { 0., 0.00627160, 0.00863282, 0.01027512, 0.01151043,
    0.01246041, 0.01318686, 0.01372531, 0.01409734, 0.01431592, 0.01438804,
    0.01431592, 0.01409734, 0.01372531, 0.01318686, 0.01246041, 0.01151043,
    0.01027512,   0.00863282,   0.00627160 };
  
  double NU9[2*nvals] = {0., -17.49285568, -16.49242250, -15.49193338,
    -14.49137675, -13.49073756, -12.48999600, -11.48912529, -10.48808848,
    -9.48683298,   8.48528137,   7.48331477,   6.48074070,   5.47722558,
    4.47213595,   3.46410162,   2.44948974,   1.41421356, 0., 0.};
  
  double MUOrg9[2*nvals] = {0.,-0.01353936,-0.01276503,-0.01199066,-0.01121624,
    -0.01044175,-0.00966718,-0.00889251,-0.00811772,-0.00734275, 0.00656756,
    0.00579204, 0.00501605, 0.00423934, 0.00346141, 0.00268119, 0.00189589,
    0.00109459, 0., 0.};
  
  double MU9[2*nvals] = { 0.,-0.03109086,-0.02931274,-0.02753453,-0.02575619,
    -0.02397771,-0.02219904,-0.02042015,-0.01864096,-0.01686139, 0.01508128,
    0.01330044, 0.01151852, 0.00973493, 0.00794853, 0.00615691, 0.00435359,
    0.00251355, 0.        , 0. };
  
  
  
};

TEST_F(ReExpConstUTest, checkA)
{
  for ( int s = -nvals + 1 ; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_a_val( nvals-1, s), A9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_a_val( nvals-1,-s), A9[-s+nvals], preclim);
  }
}

TEST_F(ReExpConstUTest, checkB)
{
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_b_val( nvals-1, s), B9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_b_val( nvals-1,-s), B9[-s+nvals], preclim);
  }
}

TEST_F(ReExpConstUTest, checkAlpha)
{
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_alpha( nvals-1, s), AL9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_alpha( nvals-1,-s), AL9[-s+nvals], preclim);
  }
  
  Constants Cst;
  set_vals(Cst.get_kappa(), 25.0);
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_alpha( nvals-1, s), AL9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_alpha( nvals-1,-s), AL9[-s+nvals], preclim);
  }
  
}

TEST_F(ReExpConstUTest, checkBeta)
{
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_beta( nvals-1, s), BEOrg9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_beta( nvals-1,-s), BEOrg9[-s+nvals], preclim);
  }
  
  Constants Cst;
  set_vals(Cst.get_kappa(), 25.0);
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_beta( nvals-1, s), BE9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_beta( nvals-1,-s), BE9[-s+nvals], preclim);
  }
}

TEST_F(ReExpConstUTest, checkNu)
{
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_nu( nvals-1, s), NU9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_nu( nvals-1,-s), NU9[-s+nvals], preclim);
  }
  
  Constants Cst;
  set_vals(Cst.get_kappa(), 25.0);
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_nu( nvals-1, s), NU9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_nu( nvals-1,-s), NU9[-s+nvals], preclim);
  }
}

TEST_F(ReExpConstUTest, checkMu)
{
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_mu( nvals-1, s), MUOrg9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_mu( nvals-1,-s), MUOrg9[-s+nvals], preclim);
  }
  
  Constants Cst;
  set_vals(Cst.get_kappa(), 25.0);
  
  for ( int s = -nvals+1; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpCo_.get_mu( nvals-1, s), MU9[ s+nvals], preclim);
    EXPECT_NEAR( ReExpCo_.get_mu( nvals-1,-s), MU9[-s+nvals], preclim);
  }
}


/*
 Class to test R and S
 */
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
  
  double RM1Real[nvals] = {-0.252453, -0.027761, 0.115623, 0.002687, 0.162679,
    -0.091828,  -0.253511,  0.141596,  0.057537,  0.118726 };
  double RM7Imag[nvals] = {-0.0,-0.103267, -0.018373, 0.230157, -0.276410,
    0.121327,  0.053113,  -0.107260,  0.066849,  -0.019357};
  
  double SN0Z[nvals] = {0.970147, 4.99775, 24.99617, 124.988519, 624.958996,
    3124.840538,  15624.347650,  78122.240044,  390613.040154, 1953072.235847};
  
  double SN0[nvals] = { 0.096100,  0.368228,  1.150332,  3.550246,  10.932428,
    33.638947,  103.470482,  318.206275,  978.479829,  3008.591425 };
  
  double SNNZ[nvals] = {0.197923, -0.599231,  1.997859,  -6.994426,  25.184186,
    -92.352278,  343.049538,  -1286.510063,  4860.364854,  0 };
  
  double SNN[nvals] = {0.368228,-10.678925,336.57133,-11139.76822,379130.294425,
    -13140242.678415,461300842.23576,-16349244619.17599,583715090846.38794, 0 };
  
  double SMMZ[nvals] = {0.19923861,  0.48877513,  0.77332760,  1.05695152,
    1.34024129, 1.62337339, 1.90641867, 2.18941109, 2.47236893, 0.00000000 };
  
  double SMM[nvals] = {1.15033201,  26.73265165,  400.10876211,  5170.65006781,
    61980.40129512,  709610.89326185,  7876260.175968,  85488801.789434,
    912350693.98094785, 0.0 };
  
  double SMLZ[nvals] = { 5.41735113,  10.98949976,  19.03910544,  29.85272945,
    43.68976929,  60.78900127,  0.0,  0.0, 0.0, 0.0 };
  
  double SML[nvals] = { 8626.68751768,  165395.13830094,  2708112.63388879,
    40130100.04533, 555039663.0819421, 7298284604.203878, 0, 0, 0, 0};
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};

TEST_F(ReExpUTest, checkR0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCo,
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
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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
  Pt testPt = Pt( 0.0, 0.0, 1.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCo,
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
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, shCon );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, bCon );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK, ReExpCoeff,
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

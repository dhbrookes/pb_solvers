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
#include <iostream>

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
  
  double shSing[10] = {0.0, 0.707106781, 1.22474487, 1.73205081, 2.23606798,
    2.73861279, 3.24037035, 3.74165739, 4.24264069, 4.74341649};
  
  double singPref0[10] = {-0.707106781,-1.22474487,-1,-1.73205081,-1.58113883,
    -1.22474487,-2.23606798,-2.12132034,-1.87082869,-1.41421356};
  
  double singPref1[10] = {0,1,0,1.58113883,1.22474487,0,2.12132034,1.87082869,
    1.41421356,0};
  
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
  
  double RM1Real[nvals] = {-0.252453478,-0.0277612693,0.115622603,0.0026871114,
    0.162679063,-0.0918281369,-0.253511354,0.14159602,0.0575370168,0.1187263};
  double RM7Imag[nvals] = {-6.51550947e-15,-0.103267311,-0.0183726894,
    0.230156648,-0.276410273,0.121326888,0.0531127252,-0.107260052,
    0.0668487462,-0.0193566816};
  
  double dRdt0real[nvals] = {-2.39498398,-0.477685728,1.04459429,-0.0529544986,
    1.4031573, 0.41740674, -2.0927581, -0.22819311, 0.45569999, -0.028397485};
  double dRdt0imag[nvals] = {0,-0.297688207,2.12864784,0.522048127,-1.8139857,
    -0.154789409, -0.428976064, -0.21665904, 1.75516662, 0.0907489077};
  
  double dRdtM1Real[nvals] = {0.5628516,-2.06104815,-0.0138369012,-0.23986527,
    0.32261620, 1.9826021, -1.0341846, -1.13418235, 0.354985375, -0.145047867};
  double dRdtM7Imag[nvals] = {-2.654661e-14,0.827169867,-1.77385735,1.44195561,
    -0.34658746,-0.152107298,-0.16301788,0.49223112,-0.39794953,0.13982327};
  
  double dRdpM1Zimag[nvals] = {0, -1, 0, 0, 0, 0, 0, 0, 0, 0};
  
  double dRdp0real[nvals] = { 0, 0.133520849, 0.106539386, -0.773127035,
    -0.18113191, 0.4777016, -0.063988443, 1.49160829, 0.34709247, -3.6890296};
  double dRdp0imag[nvals] = { 0, -0.214254385, -0.0522822197, -0.0784229506,
    -0.140109459, 1.2881752, 0.31216737, -1.5710156, -0.0901168, -1.1543848};
  
  double dRdpM1Real[nvals] = {0, -0.0173005012, 0.471225637, -0.079472085,
    -0.84123855, 0.170265856, -0.3117904, 0.94107312, 1.7728691, -3.4146876};
  double dRdpM7Imag[nvals] = {0, 0.16570801, 0.0180321104, 0.0700385421,
    -0.85523735, 1.6358568, -1.5546613, 0.79079108, -0.13884937, -0.054514482};
  
  double SN0Z[nvals] = {0.97014737,4.99774954,24.9961736,124.988519,
    624.958996,3124.84054,15624.3477,78122.24,390613.04,1953072.24};
  
  double SN0[nvals] = {0.0961000488, 0.36822761, 1.15033201, 3.55024645,
    10.9324275, 33.638947, 103.470482, 318.206275, 978.47983, 3008.59143};
  
  double SNNZ[nvals] = {0.197922989,-0.599230825,1.99785851,-6.99442583,
    25.1841864,-92.3522778,343.049538,-1286.51006,4860.36486,0 };
  
  double SNN[nvals] = {0.36822761,-10.6789251,336.571331,-11139.7682,
    379130.295,-13140242.7,461300843,-1.63492446e+10,5.83715091e+11,0};
  
  double SMMZ[nvals] = {0.19923861,  0.48877513,  0.77332760,  1.05695152,
    1.34024129, 1.62337339, 1.90641867, 2.18941109, 2.47236893, 0.00000000 };
  
  double SMM[nvals] = {1.15033201,  26.73265165,  400.10876211,  5170.65006781,
    61980.40129512,  709610.89326185,  7876260.175968,  85488801.789434,
    912350693.98094785, 0.0 };
  
  double SMLZ[nvals] = { 5.41735113,  10.98949976,  19.03910544,  29.85272945,
    43.68976929,  60.78900127,  0.0,  0.0, 0.0, 0.0 };
  
  double SML[nvals] = { 8626.68751768,  165395.13830094,  2708112.63388879,
    40130100.04533, 555039663.0819421, 7298284604.203878, 0, 0, 0, 0};
  
  
  double dSN0Z[nvals] = {-0.999549909, -9.99995464, -74.9961719, -499.977037,
    -3124.87698,-18749.3621,-109371.738,-624983.44,-3515541.28,-19530827.9};
  double dSN0[nvals] = {-0.0147291044, -0.0927621524, -0.427157004, -1.75145258,
    -6.7329656, -24.8455948, -89.130168, -313.198839, -1083.32216, -3700.72023};
  
  double dSNN[nvals] = {-0.00742097219,0.00269386658,-0.00081443389,
    0.000229928416,-6.25907958e-05,1.66587074e-05,-4.36638457e-06,
    1.13185309e-06,-2.9094711e-07,0};
  double dSMM[nvals] = {-0.00273380482,-0.000674822887,-9.03901527e-05,
    -9.60728674e-06,-9.00601584e-07,-7.79771575e-08,-6.39081501e-09,
    -5.03103237e-10,-3.84039721e-11,0};
  double dSML[nvals] = {-0.000178068655,-2.73074988e-05,-3.43353084e-06,
    -3.79874726e-07,-3.84277152e-08,-3.63796283e-09,0,0,0,0};
  double dSNNZ[nvals] = {-0.000613419148,0.000151081306,-3.01768803e-05,
    5.6494125e-06,-1.0193187e-06,1.79721615e-07,-3.11943064e-08,
    5.35330416e-09,-9.10838997e-10,0};
  double dSMMZ[nvals] = {-0.000185505058,-3.02915204e-05,-2.69622608e-06,
    -1.90297198e-07,-1.18348447e-08,-6.79347587e-10,-3.68940906e-11,
    -1.92386519e-12,-9.72507214e-14,0};
  double dSMLZ[nvals] = {-4.36830824e-06,-4.4342876e-07,-3.69119067e-08,
    -2.70365467e-09,-1.81055258e-10,-1.13459767e-11,0,0,0,0};
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};

TEST_F(ReExpUTest, checkSpecialSH)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 35.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda} );
  
  MyVector<double> singSH = ReExpTest.calc_SH_spec(1.0);
  
  for ( int s = 1; s < nvals; s++ )
    EXPECT_NEAR( singSH[s], shSing[s], preclim);
}

TEST_F(ReExpUTest, checkSingPrefactor)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 35.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda}, true);

  int ct = 0;
  for ( int s = 1; s < 5; s++ )
    for ( int m = 1; m <= s; m++ )
    {
      EXPECT_NEAR( ReExpTest.get_prefac_dR_val(s,m,0), singPref0[ct], preclim);
      EXPECT_NEAR( ReExpTest.get_prefac_dR_val(s,m,1), singPref1[ct], preclim);
      ct++;
    }
}

TEST_F(ReExpUTest, checkR0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda} );
  
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
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
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
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 5, s).real(),
                RM5Zreal[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 1, s).imag(),
                RM5Zimag, preclim);
  }
  
}

TEST_F(ReExpUTest, DISABLED_checkR)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR(ReExpTest.get_rval(nvals-1,1,s).real()/RM1Real[s],1,preclim);
    EXPECT_NEAR(ReExpTest.get_rval(nvals-1,7,s).imag()/RM7Imag[s],1,preclim);
  }
}


TEST_F(ReExpUTest, checkS0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 1.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda} );

  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_sval( 0, s, 0)/SN0Z[s],            1,preclim);
    EXPECT_NEAR( ReExpTest.get_sval( s, 0, 0)/SN0Z[s]*pow(-1.0,s),1,preclim);
  }
}


TEST_F(ReExpUTest, checkS0)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR(ReExpTest.get_sval( 0, s, 0)/SN0[s],             1, preclim);
    EXPECT_NEAR(ReExpTest.get_sval( s, 0, 0)/SN0[s]*pow(-1.0,s), 1, preclim);
  }
}

TEST_F(ReExpUTest, checkSZpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
  for ( int s = 0; s < nvals; s++ )
  {
    if ( SNNZ[s] != 0)
      EXPECT_NEAR(ReExpTest.get_sval(  s, s+1,   0)/SNNZ[s], 1, preclim);
    EXPECT_NEAR(ReExpTest.get_sval(  s, s+2,   s), SMMZ[s], preclim);
    EXPECT_NEAR(ReExpTest.get_sval(s+3, s+4, s+1), SMLZ[s], preclim);
  }
}

TEST_F(ReExpUTest, checkS)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);

  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda} );
  
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

TEST_F(ReExpUTest, checkdRdtheta0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda}, true );
  
}

TEST_F(ReExpUTest, checkdRdtheta0)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dtheta_val( nvals-1, 0, s).real(),
                dRdt0real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dtheta_val( nvals-1, 0, s).imag(),
                dRdt0imag[s], preclim);
  }
}

TEST_F(ReExpUTest, checkdRdthetaZpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 5, s).real(),
                RM5Zreal[s], preclim);
    EXPECT_NEAR( ReExpTest.get_rval( nvals-1, 1, s).imag(),
                RM5Zimag, preclim);
  }
  
}

TEST_F(ReExpUTest, checkdRdtheta)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dtheta_val( nvals-1, 1, s).real(),
                dRdtM1Real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dtheta_val( nvals-1, 7, s).imag(),
                dRdtM7Imag[s], preclim);
  }
}

TEST_F(ReExpUTest, checkdRdphi0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 0, s).real(), 0, preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 0, s).imag(), 0, preclim);
  }
  
}

TEST_F(ReExpUTest, checkdRdphi0)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 0, s).real(),
                dRdp0real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 0, s).imag(),
                dRdp0imag[s], preclim);
  }
}

TEST_F(ReExpUTest, checkdRdphiZpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 5.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 5, s).real(), 0, preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 1, s).imag(),
                dRdpM1Zimag[s], preclim);
  }
  
}

TEST_F(ReExpUTest, checkdRdphi)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );

  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 1, s).real(),
                dRdpM1Real[s], preclim);
    EXPECT_NEAR( ReExpTest.get_dr_dphi_val( nvals-1, 7, s).imag(),
                dRdpM7Imag[s], preclim);
  }
}


TEST_F(ReExpUTest, checkdSdr0Zpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 1.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCo( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCo),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR(ReExpTest.get_dsdr_val(0,s,0)/dSN0Z[s], 1.0, preclim);
    EXPECT_NEAR(ReExpTest.get_dsdr_val(s,0,0)/
                (dSN0Z[s]*pow(-1.0,s)), 1.0, preclim);
  }
}


TEST_F(ReExpUTest, checkdSdr0)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 25.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );

  for ( int s = 0; s < nvals; s++ )
  {
    EXPECT_NEAR(ReExpTest.get_dsdr_val(0,s,0)/dSN0[s],            1,preclim);
    EXPECT_NEAR(ReExpTest.get_dsdr_val(s,0,0)/dSN0[s]*pow(-1.0,s),1,preclim);
  }
  
}

TEST_F(ReExpUTest, checkdSdrZpt)
{
  Constants Cst;
  Pt testPt = Pt( 0.0, 0.0, 25.0);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 5.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    if ( dSNN[s] != 0)
      EXPECT_NEAR(ReExpTest.get_dsdr_val(  s, s+1,  0)/dSNNZ[s], 1.0, preclim);
    if ( dSNN[s] != 0)
      EXPECT_NEAR(ReExpTest.get_dsdr_val(  s, s+2,  s)/dSMMZ[s], 1.0, preclim);
    if ( dSML[s] != 0)
      EXPECT_NEAR(ReExpTest.get_dsdr_val(s+3, s+4,s+1)/dSMLZ[s], 1.0, preclim);
  }
}

TEST_F(ReExpUTest, checkdSdr)
{
  Constants Cst;
  Pt testPt = Pt( 6.9,-4.3,-0.2);
  SHCalcConstants shCon( 2*nvals );
  SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(shCon) );
  shCalc.calc_sh( testPt.theta(), testPt.phi());
  
  BesselConstants bCon( 2*nvals );
  BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(bCon) );
  
  MyMatrix<cmplx> shMat = shCalc.get_full_result();
  double kap            = Cst.get_kappa();
  double lambda         = 2.0;
  ReExpCoeffsConstants ReExpCoeff( kap, lambda, nvals);
  
  vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
  
  ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                        make_shared<ReExpCoeffsConstants> (ReExpCoeff),
                        {kap,kap}, {lambda}, true );
  
  for ( int s = 0; s < nvals; s++ )
  {
    if ( dSNN[s] != 0)
      EXPECT_NEAR( ReExpTest.get_dsdr_val(  s, s+1,  0)/dSNN[s], 1.0, preclim);
    if ( dSNN[s] != 0)
      EXPECT_NEAR( ReExpTest.get_dsdr_val(  s, s+2,  s)/dSMM[s], 1.0, preclim);
    if ( dSML[s] != 0)
      EXPECT_NEAR( ReExpTest.get_dsdr_val(s+3, s+4,s+1)/dSML[s], 1.0, preclim);
  }
}


#endif /* ReExpCalcUnitTest_h */

//
//  ElectrostaticsUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ElectrostaticsUnitTest_h
#define ElectrostaticsUnitTest_h

#include <stdio.h>
#include "Electrostatics.h"

class ElecUTest : public ::testing::Test
{
public :

protected :
  
  int vals_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -8.0 ), Pt( 0, 0, 0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -7.0 ), Pt( 0, 0, 1 ) };
    double cg[2] = { -5.0, -5.0};
    double rd[2] = {  3.6,  3.6};
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 3; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      charges[1]=cg[molInd]; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
      charges[2]=cg[molInd]; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);

      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  double min[3] = {-13.6, -13.6, -21.6};
  double bin[3] = {2.4727272727, 2.4727272727, 3.2};
  double bn2[3] = {1.36, 1.36, 1.76};
  
  double potE[1331] = {-0.00641085433,-0.0073352782,-0.00827104238,-0.00913525621,-0.00982203036,-0.0102271672,-0.0102817542,-0.0099766687,-0.00936508139,-0.00854310862,-0.00761852999,-0.00698396584,-0.00808072032,-0.00921718762,-0.0102902174,-0.0111582773,-0.0116760279,-0.0117458202,-0.0113551516,-0.0105796163,-0.00955354516,-0.00842290863,-0.00752544377,-0.00880309045,-0.0101585858,-0.0114677987,-0.0125464354,-0.0131966815,-0.0132840384,-0.0127921895,-0.0118254708,-0.0105666497,-0.00920852495,-0.00799814959,-0.00944937551,-0.0110232389,-0.0125764454,-0.0138782667,-0.0146703239,-0.0147758512,-0.0141749681,-0.0130048477,-0.0115036249,-0.00991617302,-0.00836091326,-0.00995618108,-0.0117175174,-0.0134868992,-0.0149906783,-0.0159115978,-0.0160327082,-0.0153316202,-0.0139765822,-0.0122590871,-0.0104731229,-0.00857576816,-0.0102613291,-0.0121431623,-0.0140545186,-0.0156923187,-0.016698182,-0.01682844,-0.0160595609,-0.0145806365,-0.0127204956,-0.0108068698,-0.00861718512,-0.0103211022,-0.0122276608,-0.0141677848,-0.015831399,-0.0168520188,-0.0169822767,-0.0161986412,-0.0146939027,-0.0128049942,-0.0108666429,-0.00847963365,-0.0101256903,-0.0119542506,-0.0138006141,-0.0153729193,-0.0163331236,-0.016454234,-0.0157138613,-0.0142902971,-0.0124958203,-0.0106426321,-0.00817948098,-0.00970330583,-0.0113703181,-0.0130272778,-0.0144201262,-0.0152645071,-0.0153700344,-0.0147168275,-0.0134556802,-0.0118507042,-0.0101701033,-0.00775004044,-0.00910987332,-0.0105667967,-0.0119852111,-0.0131578676,-0.0138621233,-0.0139494801,-0.0134036217,-0.0123428832,-0.0109748606,-0.00951530783,-0.00723250225,-0.00841107625,-0.00964442647,-0.010818123,-0.0117710412,-0.012337298,-0.0124070903,-0.0119679156,-0.0111075219,-0.00998078401,-0.00875326456,-0.00698396584,-0.00808072032,-0.00921718762,-0.0102902174,-0.0111582773,-0.0116760279,-0.0117458202,-0.0113551516,-0.0105796163,-0.00955354516,-0.00842290863,-0.00766564031,-0.00899305499,-0.0104103486,-0.0117878708,-0.0129286653,-0.0136185474,-0.0137113824,-0.0131897248,-0.0121668238,-0.0108403731,-0.00941737076,-0.00832147859,-0.00989997347,-0.0116391959,-0.01338335,-0.0148647093,-0.0157727557,-0.015893866,-0.0152056512,-0.013873033,-0.0121807656,-0.0104169153,-0.00890427474,-0.0107327351,-0.0128100274,-0.0149597251,-0.0168322109,-0.0179934029,-0.0181453688,-0.0172631989,-0.0155743523,-0.0134724195,-0.0113442878,-0.0093586574,-0.0114017331,-0.0137841599,-0.0163184503,-0.018574534,-0.0199843261,-0.0201636878,-0.0190898037,-0.0170509514,-0.0145559973,-0.0120943037,-0.00963111678,-0.0118125034,-0.0143995,-0.0172015078,-0.0197296869,-0.0213129933,-0.0215086237,-0.0202979762,-0.018008881,-0.0152386666,-0.0125528568,-0.00968439897,-0.011894878,-0.0145260396,-0.0173857459,-0.0199686366,-0.0215815893,-0.0217772197,-0.0205369259,-0.0181931192,-0.0153652062,-0.0126352314,-0.00951016471,-0.0116322574,-0.014131258,-0.016813507,-0.0192079953,-0.0206941659,-0.0208735276,-0.0197232651,-0.017546008,-0.0149030955,-0.012324828,-0.00913234819,-0.0110701239,-0.0133009814,-0.0156365271,-0.0176789928,-0.0189361439,-0.0190881098,-0.0181099808,-0.0162511543,-0.0139633735,-0.0116816766,-0.00859877845,-0.0102961213,-0.0121928656,-0.014117662,-0.0157599212,-0.0167601398,-0.0168812502,-0.0161008632,-0.014607345,-0.0127344353,-0.0108130631,-0.00796637615,-0.00940721196,-0.0109663024,-0.0124982228,-0.013772796,-0.0145395623,-0.0146323973,-0.0140338555,-0.0128771759,-0.0113963269,-0.00983152773,-0.00752544377,-0.00880309045,-0.0101585858,-0.0114677987,-0.0125464354,-0.0131966815,-0.0132840384,-0.0127921895,-0.0118254708,-0.0105666497,-0.00920852495,-0.00832147859,-0.00989997347,-0.0116391959,-0.01338335,-0.0148647093,-0.0157727557,-0.015893866,-0.0152056512,-0.013873033,-0.0121807656,-0.0104169153,-0.00910172469,-0.0110206547,-0.0132247285,-0.0155324724,-0.0175625621,-0.0188274797,-0.0189923417,-0.0180325785,-0.0162015325,-0.0139381938,-0.0116705456,-0.00980819348,-0.0120801917,-0.0148034375,-0.017790895,-0.0205220193,-0.0222467247,-0.0224604794,-0.0211531968,-0.0186879676,-0.0157210922,-0.0128743812,-0.0103685343,-0.0129558327,-0.0161781579,-0.0198735067,-0.0233791519,-0.0256082397,-0.0258627259,-0.0241789216,-0.0210166292,-0.0172993924,-0.0138809903,-0.0107091283,-0.0135064599,-0.0170828624};
  
  double pot2DZ[121] = {-0.0073352782,-0.00808072032,-0.00880309045,-0.00944937551,-0.00995618108,-0.0102613291,-0.0103211022,-0.0101256903,-0.00970330583,-0.00910987332,-0.00841107625,-0.00808072032,-0.00899305499,-0.00989997347,-0.0107327351,-0.0114017331,-0.0118125034,-0.011894878,-0.0116322574,-0.0110701239,-0.0102961213,-0.00940721196,-0.00880309045,-0.00989997347,-0.0110206547,-0.0120801917,-0.0129558327,-0.0135064599,-0.0136199697,-0.013268151,-0.0125243334,-0.0115248874,-0.0104099505,-0.00944937551,-0.0107327351,-0.0120801917,-0.0133934632,-0.0145128563,-0.0152359676,-0.0153892865,-0.0149263359,-0.0139625054,-0.012703119,-0.0113417621,-0.00995618108,-0.0114017331,-0.0129558327,-0.0145128563,-0.0158793674,-0.0167852522,-0.0169811431,-0.0163972572,-0.0152036067,-0.0136873753,-0.012096362,-0.0102613291,-0.0118125034,-0.0135064599,-0.0152359676,-0.0167852522,-0.0178302007,-0.0180558781,-0.0173741132,-0.0160057492,-0.014305209,-0.0125579088,-0.0103211022,-0.011894878,-0.0136199697,-0.0153892865,-0.0169811431,-0.0180558781,-0.0182815555,-0.0175700041,-0.016159068,-0.0144187189,-0.0126402834,-0.0101256903,-0.0116322574,-0.013268151,-0.0149263359,-0.0163972572,-0.0173741132,-0.0175700041,-0.016915147,-0.0156170863,-0.0139996936,-0.0123268863,-0.00970330583,-0.0110701239,-0.0125243334,-0.0139625054,-0.0152036067,-0.0160057492,-0.016159068,-0.0156170863,-0.0145315476,-0.0131472607,-0.0116791509,-0.00910987332,-0.0102961213,-0.0115248874,-0.012703119,-0.0136873753,-0.014305209,-0.0144187189,-0.0139996936,-0.0131472607,-0.0120291201,-0.0108060983,-0.00841107625,-0.00940721196,-0.0104099505,-0.0113417621,-0.012096362,-0.0125579088,-0.0126402834,-0.0123268863,-0.0116791509,-0.0108060983,-0.00982136892};
  
  string dx[9] = {"# Data from PBAM Electrostat run",
    "# My runname is /Users/lfelberg/Desktop/test.dx and units internal",
    "object 1 class gridpositions counts 11 11 11",
    "origin -13.6 -13.6 -21.6",
    "delta 2.47273 0.0e+00 0.0e+00",
    "delta 0.0e00 2.47273 0.0e+00",
    "delta 0.0e00 0.0e+00 3.2",
    "object 2 class gridconnections counts 11 11 11",
    "object 3 class array type double rank 0 items 1331 data follows"};
  
  string potx[7] = {"# Data from PBAM Electrostat run",
    "# My runname is /Users/lfelberg/Desktop/pot_x_1.00.dat",
    "units internal", "grid 11 11", "axis x 1",
    "origin -13.6 -21.6", "delta 2.47273 3.2"};
  string poty[7] = {"# Data from PBAM Electrostat run",
    "# My runname is /Users/lfelberg/Desktop/pot_y_10.00.dat",
    "units internal", "grid 11 11", "axis y 10",
    "origin -13.6 -21.6", "delta 2.47273 3.2"};
  string potz[7] = {"# Data from PBAM Electrostat run",
    "# My runname is /Users/lfelberg/Desktop/pot_z_-17.23.dat","units internal",
    "grid 11 11","axis z -17.23","origin -13.6 -13.6","delta 2.47273 2.47273"};
  
  
} ; // end EnergyForceUTest

TEST_F(ElecUTest, detailsCheck)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  vector<double>  mins = EstatTest.get_mins();
  vector<double>  maxs = EstatTest.get_maxs();
  vector<int>     npts = EstatTest.get_npts();
  vector<double>  bins = EstatTest.get_bins();
  
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR( mins[i]/min[i], 1, preclim);
    EXPECT_NEAR(   maxs[i]/13.6, 1, preclim);
    EXPECT_NEAR(   npts[i]/11.0, 1, preclim);
    EXPECT_NEAR( bins[i]/bin[i], 1, preclim);
  }
}

TEST_F(ElecUTest, detailsCheck2)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 20);
  vector<double>  mins = EstatTest.get_mins();
  vector<double>  maxs = EstatTest.get_maxs();
  vector<int>     npts = EstatTest.get_npts();
  vector<double>  bins = EstatTest.get_bins();
  
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR( mins[i]/min[i], 1, preclim);
    EXPECT_NEAR(   maxs[i]/13.6, 1, preclim);
    EXPECT_NEAR(   npts[i]/20.0, 1, preclim);
    EXPECT_NEAR( bins[i]/bn2[i], 1, preclim);
  }
}

TEST_F(ElecUTest, checkPot)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  vector<vector<vector<double > > > esp = EstatTest.get_potential();
  
  int ct = 0;
  for (int i = 0; i < esp.size(); i++)
    for (int j = 0; j < esp[0].size(); j++)
      for (int k = 0; k < esp[0][0].size(); k++)
      {
        if (ct < 300)
        {
          EXPECT_NEAR( esp[i][j][k]/potE[ct], 1, preclim);

        } else
          return;
        ct++;
      }
}

TEST_F(ElecUTest, printDX)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  EstatTest.print_dx("/Users/lfelberg/Desktop/test.dx");
  
  string inputLine;
  ifstream fin("/Users/lfelberg/Desktop/test.dx");
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 9) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == dx[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, checkPOTZ)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = -17.23;
  char pot[50];
  sprintf(pot, "/Users/lfelberg/Desktop/pot_z_%.2f.dat", val);
  EstatTest.print_grid(Zdim, val, string(pot));
  
  vector<vector<double > > p2d = EstatTest.get_pot2d();
  
  int ct = 0;
  for (int i = 0; i < p2d.size(); i++)
    for (int j = 0; j < p2d[0].size(); j++)
    {
      if (ct < 121) EXPECT_NEAR( p2d[i][j]/pot2DZ[ct], 1, preclim);
      else          return;
      ct++;
    }
}

TEST_F(ElecUTest, printPOTX)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = 1;
  char pot[50];
  sprintf(pot, "/Users/lfelberg/Desktop/pot_x_%.2f.dat", val);
  EstatTest.print_grid(Xdim, val, string(pot));
  
  string inputLine;
  ifstream fin(pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == potx[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOTY)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = 10;
  char pot[50];
  sprintf(pot, "/Users/lfelberg/Desktop/pot_y_%.2f.dat", val);
  EstatTest.print_grid(Ydim, val, string(pot));
  
  string inputLine;
  ifstream fin(pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == poty[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOTZ)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
  Electrostatic EstatTest( ASolvTest, 11);
  
  double val = -17.23;
  char pot[50];
  sprintf(pot, "/Users/lfelberg/Desktop/pot_z_%.2f.dat", val);
  EstatTest.print_grid(Zdim, val, string(pot));
  
  string inputLine;
  ifstream fin(pot);
  getline(fin,inputLine);
  
  int ct = 0;
  while ( (ct < 7) && (!fin.eof()))
  {
    EXPECT_TRUE( inputLine == potz[ct]);
    getline(fin,inputLine);
    ct++;
  }
}

TEST_F(ElecUTest, printPOT)
{
  const int vals = 5;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ASolvTest->solve_A(1E-12); ASolvTest->solve_gradA(1E-12);
  
//  Electrostatic EstatTest( ASolvTest, 111);
//  double val = 0; char pot[50];
//  sprintf(pot, "/Users/lfelberg/Desktop/pot_x_%.2f.dat", val);
//  EstatTest.print_grid(Xdim, val, string(pot));
//  
//  EXPECT_TRUE( 0 == 0);
}



#endif /* ElectrostaticsUnitTest_h */

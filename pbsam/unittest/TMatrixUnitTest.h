//
//  TMatrixUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 6/16/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef TMatrixUnitTest_h
#define TMatrixUnitTest_h

#include "Solver.h"

/*
 Class for unit testing transform class
 */
class TMatrixUTest : public ::testing::Test
{
  protected :
  
  virtual void SetUp()   {}
  virtual void TearDown() {}
  
  vector<vector<vector<double> > > localXSphre = {{{0.0016035307,0.00015733459,-0.000220676033,-0.000277486764,-5.84156321e-05,-0.000281432114,-4.15405527e-05,3.62790353e-05,-5.48851976e-06,0.000158864088,7.28883296e-05,1.64229817e-05,7.03977617e-05,3.37574114e-05,2.47195018e-05,},{0.00971106004,0.00179553222,0.00231271522,-0.000519470092,0.000739615662,0.000315991496,-0.000359412724,-4.84715361e-05,0.000160339366,-5.99437478e-05,-2.07536799e-05,-0.000104825193,1.79579618e-05,-9.78039148e-06,-4.17535012e-05,},{0.0151132106,0.00461152795,0.00289852255,0.000461290198,0.0013728625,0.000207432747,-0.000291219344,0.000291302654,5.37355373e-05,-0.000197095971,-0.000114813533,-1.01429672e-05,-2.5730573e-05,-0.000157837744,-0.000121378356,},{0.00029268089,0.000140412729,-2.39183684e-05,4.29173001e-05,-2.09584277e-05,-2.27743778e-05,-2.23233215e-06,-1.41658243e-05,-2.26910752e-05,3.57985293e-06,-1.37795899e-05,-7.96378244e-06,-1.36293037e-05,4.07150193e-06,2.54720501e-06,},{-0.000381219794,-0.000232917003,3.73470718e-05,-0.0001366731,3.81198649e-05,-1.2437576e-05,-7.62653811e-05,3.05187848e-05,-1.40307701e-05,6.1542867e-07,-4.0162007e-05,2.13248058e-05,-1.19537416e-05,7.38483322e-07,5.46163835e-07,},{0.00313899898,0.000218642586,-0.000671110512,-0.000340566087,-3.89955752e-05,-6.05600457e-05,-5.18905088e-05,9.11694164e-05,-4.74906202e-05,0.000136478965,4.29683315e-05,3.61287285e-07,-8.67547995e-06,1.06813077e-05,-7.81472483e-05,},{0.000868708013,0.000235136379,0.000288156023,-3.72354702e-05,0.000109963709,0.000124045136,-3.34530568e-05,-1.52031815e-05,5.11125273e-05,5.54175769e-05,1.75274282e-05,-1.70523381e-05,-9.89459086e-06,2.36246979e-05,2.29331389e-05,},{-0.000861703097,-0.000393875261,2.93251481e-05,-0.000142508135,1.07016375e-06,2.2238505e-05,-2.65218696e-05,-8.71782516e-06,2.75867055e-05,-1.50419614e-05,8.69284952e-06,-6.08864449e-06,1.96467325e-05,-7.84397232e-06,-1.83047419e-06,},{0.000239480494,-0.000153969001,3.07033843e-05,8.56694021e-05,-2.53435944e-05,-3.35189055e-06,-4.33247026e-05,1.23528816e-05,6.37954393e-06,-4.60440456e-06,2.19321777e-05,-1.82703758e-06,-9.22134008e-06,3.71751294e-06,-2.83412893e-06,},{-2.99786257e-07,0.000105896716,0.000120047161,6.98322605e-06,-1.83971204e-05,0.000190955935,-7.69220958e-05,-1.57781237e-05,-0.000134456438,3.43811413e-05,4.46501404e-05,7.03769887e-06,3.4340659e-05,-3.32536746e-05,-6.78471848e-05,},{0.000136980927,-7.10262562e-05,8.0956837e-05,-5.99206609e-05,-2.99985153e-05,-1.27010197e-05,2.79005685e-05,-1.60856492e-05,2.72345647e-05,-6.45861201e-05,3.52483067e-05,1.20321371e-05,1.8326764e-05,4.20647189e-05,-2.15782756e-05,},{-0.00432924649,-0.000838186199,0.00244666673,0.00142217286,0.000855321281,-0.000864637893,0.000864888668,-0.000745594496,-0.000527093356,-4.31095854e-06,-0.000408961199,-0.000600549758,0.000142633823,0.000176030361,0.000327383278,},{},{-0.00252080082,0.00158346761,0.000945977732,-0.000261380204,-0.00102191319,2.10100845e-05,-0.000744016141,0.000564182936,-3.50487374e-05,-0.000441365093,0.000992266933,8.08802341e-05,3.79454029e-05,0.000727604526,0.000420303895,},{},{-0.00647012024,0.00544427554,0.00191256524,-0.00393465145,-0.0027509801,-0.000669716952,0.00225871806,0.0029732723,0.00122343933,0.000241223281,-0.000731555605,-0.00265201405,-0.00162548283,-0.000511111066,-8.79423972e-05,}}, {{-0.00278782486,0.000679110611,0.000363809117,-7.31794743e-05,-0.000153231653,8.32567182e-06,-2.10152335e-05,3.96141313e-05,-8.51122633e-07,-2.4033008e-05,1.2802049e-05,-5.93645341e-06,-2.0292142e-06,1.31266764e-05,8.10243342e-06,},{0.00196042605,-0.00077443132,-0.000605950613,7.64769137e-05,0.000341397206,0.000177897254,9.34367839e-05,-7.21963428e-05,-7.28902502e-05,-4.63575195e-05,-4.99764137e-05,-1.28820629e-05,-1.72322494e-05,-4.58170443e-06,1.76330731e-05,},{0.0171977281,-0.00235203871,0.00770292243,-0.00310494141,-0.00116521722,0.00361887749,0.000450561305,-0.00158040959,-0.000237570564,0.00150451943,0.000616610215,1.38000937e-05,-0.00109565625,5.75192096e-05,0.000451496832,},{0.0136195272,0.000798362336,0.00124938119,-0.000835897003,0.000129556883,-0.000804326076,-0.00019178366,-8.15747712e-05,-0.000137801618,-0.000296290624,4.37681966e-05,-2.7117341e-05,3.74182307e-05,-6.20576793e-05,5.93273083e-06,},{0.000201577368,3.89155627e-06,-2.01378303e-05,-1.48077773e-05,-1.37201669e-06,-1.22641883e-05,-1.06716987e-06,1.85512182e-06,-2.32391459e-07,5.3161274e-06,1.60508272e-06,3.21225204e-07,1.12498757e-06,4.82472372e-07,-1.93795337e-07,},{-0.000253893197,0.000121425006,-3.40287887e-05,0.000138518724,-0.000118060733,6.86982511e-05,-1.54706899e-05,-8.6246284e-05,4.21985298e-05,-1.31527764e-05,-0.000107415037,5.02605077e-06,-7.66352423e-06,2.85726447e-05,-6.19749047e-06,},{0.00186332079,-0.000119204966,-0.000161496462,-5.28188803e-05,1.91190273e-05,-3.9264e-05,1.05401962e-05,4.44548284e-06,4.47588799e-06,1.66905383e-05,1.52855341e-06,-1.42633255e-06,1.24260654e-06,-2.67940433e-06,-1.7071779e-06,},{0.000925995562,2.73952976e-05,0.000199614234,-0.000105813077,1.31035176e-05,-4.2225979e-05,-2.94347728e-05,-2.5164287e-05,-1.5274777e-05,-6.57762723e-05,8.12060369e-06,-1.08081786e-05,-5.19735244e-07,-1.97551939e-05,-3.23347276e-05,},{0.000120849533,-3.5201549e-05,-4.88793165e-07,8.88056023e-06,4.24410844e-07,-1.65492553e-06,-1.90057476e-06,-2.48028118e-07,9.60957509e-07,-2.46531059e-09,3.05557015e-07,1.13963593e-07,-4.17188281e-07,-1.83865246e-08,3.02971889e-08,},{9.22344579e-05,4.96572121e-05,0.000114296042,-7.66706915e-05,-7.40329697e-05,1.9293386e-05,5.81381729e-05,3.55239063e-05,-1.54582065e-05,3.64486702e-06,-3.32144626e-05,-1.44546425e-05,8.72380395e-06,-3.5696235e-06,5.33249208e-07,},{5.91419862e-05,-1.72889059e-05,5.91091306e-06,2.28624824e-06,-1.19693609e-06,-9.56371081e-07,1.33229863e-06,-8.51262284e-08,1.672175e-06,-1.24450076e-06,-6.81068533e-07,8.44326781e-08,-6.24245308e-07,3.40410831e-07,-9.30469662e-08,},{-0.00278393419,0.000753290972,0.000836577066,9.07774134e-05,-0.000368599042,-0.000249509539,-0.000173882012,2.27613364e-05,0.000121609588,4.36921241e-05,4.32645894e-05,5.83829071e-05,-2.05989583e-05,-1.80502873e-05,8.43550564e-06,},{},{-0.00119371096,0.000406741593,0.000159321112,-0.000106264325,-9.39411967e-05,-1.24860074e-05,1.42508931e-05,3.99389448e-05,9.5068343e-06,-2.53328939e-06,},{-0.00265108177,0.000729866989,0.000225051472,-0.000145205584,-0.000107281643,2.1443137e-05,9.51078483e-06,3.59658783e-05,-1.30198314e-05,-1.21456204e-05,}}};
  
  vector<vector<vector<vector<double> > > > Hin = {{{{-0.00585,-0.0111658741,0.0508634184,-0.052510966,-0.00503811243,-0.0326385457,0.00908195863,-0.0165499472,0.00417366401,-0.043845211,0.020661115,0.00372575487,0.0110565917,0.00663396028,-0.0150460049}, {0.0,0.0,0.0882629907,0.0,-0.00874260686,0.0563207937,0.0,-0.0287190261,-0.00720203871,-0.000106503325,0.0,0.00646528051,-0.0190791595,1.61143899e-05,-0.026256467}}, {{-0.2011,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {0.0,0.0,0,0.0,0,0,0.0,0,0,0,0.0,0,0,0,0}}}, {{{0.00438225415,0.000184416296,0.00285261912,-0.000896306288,0.00109196276,0.000558332113,5.89605634e-05,-7.61231311e-05,0.00119903431,6.96383244e-05,-0.0004028978,0.000115874699,-0.000167683669,0.00136333431,     0.000640611311},{0.0,0.0,-0.000680999907,0.0,0.00124067469,-0.00103087442,0.0,0.000441185826,0.000863508264,0.000234523443,0.0,-0.000365710005,0.000678223951,0.000253506416,9.38494831e-05}}}};
  
  vector<vector<vector<double> > > Hout_re = {{{0.000782700702,-0.000791370539,-0.000240787154,0.000315397505,0.000271212625,0.00011115778,-3.46729293e-05,-0.000140086331,-8.64384885e-05,-3.8480857e-05,-3.27199469e-05,3.88124732e-05,4.20893571e-05,2.30189355e-05,9.23807324e-06,},{-0.0126186235,0.00429416833,0.00249780646,-0.000663089408,-0.00112838681,-0.000368626104,-6.40250261e-05,0.000324370933,0.000210481912,3.40662188e-05,8.71256387e-05,-5.19089875e-05,-8.00924331e-05,-2.29175302e-05,3.34503709e-06,}},{{0.000219711312,2.4359843e-05,1.56976434e-05,3.9595817e-07,2.39467155e-06,-9.55848593e-08,-2.49210848e-07,2.15893926e-07,4.31303632e-08,-1.31529139e-07,-4.87608466e-08,4.58566939e-09,2.37357946e-08,-2.21373536e-08,-1.22779539e-08,}}};
  vector<vector<vector<double> > > Hout_im = {{{0,0,0.000128274702,0,5.0240829e-05,-2.13656528e-05,0,-5.82794396e-05,-4.13222261e-05,-1.07636756e-05,0,2.25209434e-05,3.99603276e-05,2.35503497e-05,9.11933724e-06,},{0,0,0.00113296594,0,-0.000511818607,-0.000421027815,0,0.000147129582,0.000240402779,0.00011279649,0,-2.35451048e-05,-9.14779007e-05,-7.58821218e-05,-2.50925929e-05,}},{{0,0,-1.4608819e-05,0,-1.8946681e-06,-1.86630436e-06,0,-7.42577582e-08,-3.4554469e-07,-9.07496691e-08,0,2.28743061e-08,-3.94441672e-08,-2.99867886e-08,3.41056879e-09,}}};
  
};

TEST_F(TMatrixUTest, xfor_numeric_test)
{
  int pol = 5;
  double kap = 0.0;
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<MoleculeSAM> > mols;
  mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                     pqr.get_atom_pts(), pqr.get_radii(),
                                     pqr.get_cg_centers(), pqr.get_cg_radii()));
  auto sys = make_shared<System>(mols);
  
  auto cst = make_shared<Constants> ();
  cst->set_salt_concentration(0.0);
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  IEMatrix ieMatTest(0, sys->get_MoleculeSAM(0), SHCalcTest, pol,
                     _expcons, true, 0, true);

  auto ReExp = make_shared<ReExpCoeffsConstants>(kap,sys->get_lambda(),pol);
  
  // Analytical H matrix initialized
  auto hmat = make_shared<HMatrix>(0, sys->get_Ns_i(0), pol, kap);
  hmat->init(sys->get_MoleculeSAM(0), SHCalcTest, 4.0);
  
  // Local H, numeric
  LHMatrix lhmt(0, sys->get_Ns_i(0), pol, kap);
  lhmt.init(sys->get_MoleculeSAM(0), hmat, SHCalcTest, BesselCal, _expcons);
  
  // Transform class
  TMatrix tmat( pol, sys, SHCalcTest, cst, BesselCal, ReExp);
  
  vector<int> mySphs = {0, 8};
  
  for (int i = 0; i < mySphs.size(); i++)
  {
    int sphct = 0;
    for (int j = 0; j < sys->get_Ns_i(0); j++)
    {
      if ( j == mySphs[i] ) continue;
      if ( localXSphre[i][sphct].size() == 0 ) continue;
      int expct = 0;
      MyMatrix<cmplx> out = tmat.re_expandX_numeric(lhmt.get_mat(), 0,
                                                    mySphs[i], 0, j, kap);
      for (int n = 0; n < pol; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          EXPECT_NEAR(localXSphre[i][sphct][expct],out(n, m+pol).real(),
                      preclim);
          expct++;
        }
      }
      sphct++;
    }
  }
}

TEST_F(TMatrixUTest, xforIntra_analytic_test)
{
  int pol = 5;
  double kap = 0.21053961;
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<MoleculeSAM> > mols;
  mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                     pqr.get_atom_pts(), pqr.get_radii(),
                                     pqr.get_cg_centers(), pqr.get_cg_radii()));
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> ();
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  IEMatrix ieMatTest(0, sys->get_MoleculeSAM(0), SHCalcTest, pol,
                     _expcons, true, 0, true);
  
  auto ReExp = make_shared<ReExpCoeffsConstants> (kap, sys->get_lambda(), pol);
  
  // Transform class
  TMatrix tmat( pol, sys, SHCalcTest, cst, BesselCal, ReExp);
  
  vector<int> mySphs = {3, 16};
  vector<vector<int> > myXF = {{14, 16}, {8}};
  
  for (int i = 0; i < mySphs.size(); i++)
  {
    int sphct = 0;
    for (int j = 0; j < myXF[i].size(); j++)
    {
      // Analytical H matrix initialized
      int ct = 0;
      MyMatrix<cmplx> hin(pol, 2*pol+1);
      for (int n = 0; n < pol; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          hin(n, m+pol) = complex<double> (Hin[i][j][0][ct], Hin[i][j][1][ct]);
          if (m > 0) hin(n, -m+pol) = conj(hin(n, m+pol));
          ct++;
        }
      }
      
      int expct = 0;
      MyMatrix<cmplx> out = tmat.re_expandX(hin, 0,
                                            mySphs[i], 0, myXF[i][j]);
      for (int n = 0; n < pol; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          EXPECT_NEAR(Hout_re[i][sphct][expct],out(n, m+pol).real(),
                      preclim);
          EXPECT_NEAR(Hout_im[i][sphct][expct],out(n, m+pol).imag(),
                      preclim);
          expct++;
        }
      }
      sphct++;
    }
  }
}

#endif /* TMatrixUnitTest_h */

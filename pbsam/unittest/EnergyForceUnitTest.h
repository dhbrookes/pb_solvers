//
//  EnergyForceUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 7/25/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef EnergyForceUnitTest_h
#define EnergyForceUnitTest_h

#include "PhysCalcSAM.h"
#include "Solver.h"

/*
 Class for unit testing energy calculations
 */
class EnergyUTest : public ::testing::Test
{
protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}
  
  vector<double> en3 = {-0.000883187183, -0.000457557568,-0.000425629615};
  
};

TEST_F(EnergyUTest, two_mol_test)
{
  int pol(3), nmol(2);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<BaseMolecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
  auto sys = make_shared<SystemSAM>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_moli(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "grad_test/mpol.";
  vector<vector<string> > imat_loc(sys->get_n());
  vector<vector<vector<string > > > exp_loc(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      exp_loc[i][k].resize(2);
      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
      exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+ ".H.exp";
      exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+ ".F.exp";
    }
  }
  
  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                  true, true, imat_loc, exp_loc);
  solvTest.update_LHN_all();
  
  EnergyCalcSAM ecal (sys->get_n());
  shared_ptr<vector<double> > ene;
  ecal.calc_all_energy(solvTest.get_all_H(), solvTest.get_all_LHN());
  ene = ecal.get_omega();
  for (int i = 0; i < sys->get_n(); i++)
  {
    EXPECT_NEAR((*ene)[i]/-0.00042138573, 1.0, preclim);
  }
}


TEST_F(EnergyUTest, three_mol_test)
{
  int pol(3), nmol(3);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<BaseMolecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
  mols[2]->translate(Pt(-22.28786458,-14.35779167,-15.15628125), 1e14);
  auto sys = make_shared<SystemSAM>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_moli(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "grad_test/mpol3.";
  vector<vector<string> > imat_loc(sys->get_n());
  vector<vector<vector<string > > > exp_loc(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      exp_loc[i][k].resize(2);
      imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
      exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+".H.exp";
      exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+".F.exp";
    }
  }
  
  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                  true, true, imat_loc, exp_loc);
  solvTest.update_LHN_all();
  
  EnergyCalcSAM ecal (sys->get_n());
  shared_ptr<vector<double> > ene;
  ecal.calc_all_energy(solvTest.get_all_H(), solvTest.get_all_LHN());
  ene = ecal.get_omega();
  for (int i = 0; i < sys->get_n(); i++)
  {
//    cout << "This is energy "<< setprecision(9) << ene[i] << endl;
    EXPECT_NEAR((*ene)[i]/en3[i], 1.0, preclim);
  }
}


/*
 Class for unit testing solver matrices
 */
class ForceUTest : public ::testing::Test
{
  protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}

  vector<vector<vector<double> > > forIk2 = {{{-0.000444176673,-0.00029853237,-0.000530048752},{-5.61536443e-05,-5.46468477e-05,-4.89727287e-05},{0.000147779405,0.00114560544,0.000728103886},{-8.62847785e-05,0.000568286339,-0.000474743491},{0.000163021734,2.89614299e-05,9.0441767e-05},{0.000170631725,0.000215503896,0.000120660003},{-8.77086934e-05,-1.9920445e-05,-9.54192069e-05},{0.000383563835,0.000752832542,0.00061458346},{4.96561803e-05,0.000280000544,0.000118658649},{-0.000104386801,-8.80826184e-05,-0.000179106285},{-2.17720452e-05,2.59176133e-05,1.8181021e-05},{5.77513438e-05,5.95953743e-06,6.08268383e-05},{-6.64198405e-06,-4.17539093e-05,1.25891743e-06},{0,0,0},{-3.74833078e-05,-3.17878601e-05,-4.90512409e-05},{0,0,0},{-0.000114564351,-7.42639163e-05,-0.00016732068}},{{-0.000942837748,-0.00151625162,-0.00147534227},{-0.000112532561,-0.000108727631,-0.000182055666},{0.000345766375,0.000126536859,0.000397627506},{6.05343637e-05,9.74518861e-05,8.91149859e-05},{1.29014838e-05,3.40606516e-05,4.28721667e-05},{1.83338313e-05,1.54310128e-05,2.11771342e-05},{-6.49940585e-05,-0.000147060415,-0.000113951157},{4.39441985e-05,3.89276746e-05,5.86523514e-05},{6.97091269e-05,1.59124907e-05,0.000100899852},{-0.000247856795,-0.000421276713,-0.000283140494},{-2.62320988e-05,-3.96521948e-05,-8.99705811e-05},{-9.83514402e-06,-2.56603608e-05,-8.82923616e-06},{-0.000104116334,-0.00020462616,-0.000336569414},{0,0,0},{-0.000151615085,-0.000222839559,-0.000202056249},{0,0,0},{-0.000285200902,-0.000574989672,-0.000331577164}}};
  vector<vector<double> > for2 = {{1.32319456e-05,0.00241407937,0.000208052158},
    {-0.00139403135, -0.00293276375, -0.00231314824}};
  vector<vector<double> > tor2 = {{-0.0103130859,0.000723565671,0.00870034099},
    {-0.0140884639, 0.00121824127, 0.0079666432}};
  
  vector<vector<vector<double> > > forIk3 = {{{-0.00142869113,-0.00184783503,-0.0020568363},{-0.000185865177,-0.000182622931,-0.000253792692},{0.000532392856,0.00131047448,0.00119173142},{1.03220592e-05,0.000735576936,-0.000341512161},{0.000138151359,4.67005861e-05,0.000103267929},{0.000245917556,0.000268857088,0.000208497434},{-0.000181634419,-0.00017140386,-0.00023646897},{0.000462821813,0.000843108838,0.000734394356},{9.59575926e-05,0.000256026756,0.00019518453},{-0.000373016628,-0.000527868019,-0.000495993106},{-7.38274658e-05,-6.53744519e-05,-0.00014066043},{1.63023283e-05,-3.17621523e-05,1.75500746e-05},{-6.67039865e-05,-0.000193323607,-0.000251845122},{0,0,0},{-0.000192374143,-0.000256650248,-0.000255274952},{0,0,0},{-0.000407818354,-0.000654170991,-0.000509012792},},{{-0.000998241253,-0.00164418264,-0.00156529623},{-0.000120040617,-0.000118544596,-0.000192296127},{0.000379581882,0.000148304592,0.00043499646},{6.99867358e-05,0.000112327938,0.000102933731},{1.39176707e-05,3.86363067e-05,4.64922817e-05},{1.98367854e-05,1.73558909e-05,2.29778474e-05},{-6.84766665e-05,-0.000161587141,-0.000121170152},{5.03183052e-05,4.56294881e-05,6.71844489e-05},{7.69503414e-05,1.98574341e-05,0.000110372748},{-0.000263130867,-0.000456345696,-0.000301988213},{-2.56389994e-05,-4.25470925e-05,-9.49148182e-05},{-1.02943024e-05,-2.82467916e-05,-9.39045363e-06},{-0.000112111044,-0.000225469403,-0.000359703184},{0,0,0},{-0.000160856479,-0.000239768033,-0.000214590699},{0,0,0},{-0.000300561649,-0.00062332628,-0.000351072179},},{{-0.000509779661,-0.000344074656,-0.00060777769},{-6.36894738e-05,-6.06794015e-05,-5.72128641e-05},{0.000209214548,0.0012467362,0.000832171251},{-4.09310153e-05,0.000611579341,-0.000424977643},{0.000177185597,3.52008327e-05,0.000102370852},{0.000186930981,0.000230165406,0.000136831806},{-9.82889907e-05,-2.56648262e-05,-0.000108161702},{0.000421652073,0.000794153822,0.000664244404},{6.40515324e-05,0.00030343604,0.000139056884},{-0.000120793092,-0.000100332576,-0.000201601907},{-2.24065315e-05,2.76834762e-05,2.04709744e-05},{6.11028304e-05,7.43830813e-06,6.40645367e-05},{-1.50126427e-05,-4.88306128e-05,-4.80829345e-06},{0,0,0},{-4.34941229e-05,-3.65063989e-05,-5.66146484e-05},{0,0,0},{-0.000133137174,-8.69689886e-05,-0.000192143802}}};

  vector<vector<double> > for3 = {{-0.00140806574, -0.000470266608,
    -0.00209077078},{-0.00144876016, -0.00315790603, -0.00242546454},
    {7.26048576e-05, 0.00255333596, 0.000305912159}};
  vector<vector<double> > tor3 = {{-0.0251299755, 0.00202261569, 0.0173416251},
    {-0.0152504508,0.00121797971,0.00867183059},{-0.0111485196,0.000571963184,
      0.00947117611}};
  
};

TEST_F(ForceUTest, two_mol_test)
{
  int pol(3), nmol(2);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<BaseMolecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
  auto sys = make_shared<SystemSAM>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_moli(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "grad_test/mpol.";
//  string estart = test_dir_loc + "spol_test/test_0.00_p3.0.";
  vector<vector<string> > imat_loc(sys->get_n());
  vector<vector<vector<string > > > exp_loc(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      exp_loc[i][k].resize(2);
      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
      exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+ ".H.exp";
      exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+ ".F.exp";
//      exp_loc[i][k][0] = estart+to_string(k)+ ".H.exp";
//      exp_loc[i][k][1] = estart+to_string(k)+ ".F.exp";
    }
  }
  
  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                  true, true, imat_loc, exp_loc);
  solvTest.update_LHN_all();
  
  GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
                       solvTest.get_all_F(), solvTest.get_all_H(),
                       solvTest.get_IE(), solvTest.get_interpol_list(),
                       solvTest.get_precalc_sh(), _expcons, pol);
  gsolvTest.solve(1e-16, 100);
  
  auto focal = make_shared<ForceCalcSAM> (nmol, sys->get_all_Ik(),
                                       cst->get_dielectric_water(),
                                       SHCalcTest, BesselCal);
  focal->calc_all_f(solvTest.get_all_H(), solvTest.get_all_LHN(),
                    gsolvTest.get_gradH_all(),gsolvTest.get_gradLHN_all());
  shared_ptr<vector<Pt> > fo = focal->get_F();
  
  TorqueCalcSAM tocal(nmol);
  tocal.calc_all_tau(sys, focal);
  shared_ptr<vector<Pt> > to = tocal.get_tau();
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    vector<Pt> fI = focal->get_all_fIk(i);
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      for (int d = 0; d < 3; d++)
      {
        if (k==0)
        {
          EXPECT_NEAR(for2[i][d]/(*fo)[i].get_cart(d), 1.0, preclim);
          EXPECT_NEAR(tor2[i][d]/(*to)[i].get_cart(d), 1.0, preclim);
        }
        if (fabs(forIk3[i][k][d]) > 1e-11)
          EXPECT_NEAR(forIk2[i][k][d]/fI[k].get_cart(d), 1.0, preclim);
      }
    }
  }
}


TEST_F(ForceUTest, three_mol_test)
{
  int pol(3), nmol(3);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<BaseMolecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
  mols[2]->translate(Pt(-22.28786458,-14.35779167,-15.15628125), 1e14);
  auto sys = make_shared<SystemSAM>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_moli(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "grad_test/mpol3.";
  vector<vector<string> > imat_loc(sys->get_n());
  vector<vector<vector<string > > > exp_loc(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      exp_loc[i][k].resize(2);
      imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
      exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+ ".H.exp";
      exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+ ".F.exp";
    }
  }
  
  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                  true, true, imat_loc, exp_loc);
  solvTest.update_LHN_all();

  GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
                       solvTest.get_all_F(), solvTest.get_all_H(),
                       solvTest.get_IE(), solvTest.get_interpol_list(),
                       solvTest.get_precalc_sh(), _expcons, pol);
  gsolvTest.solve(1e-16, 85);

  auto focal = make_shared<ForceCalcSAM> (nmol, sys->get_all_Ik(),
                                       cst->get_dielectric_water(),
                                       SHCalcTest, BesselCal);
  focal->calc_all_f(solvTest.get_all_H(), solvTest.get_all_LHN(),
                   gsolvTest.get_gradH_all(),gsolvTest.get_gradLHN_all());
  shared_ptr<vector<Pt> > fo = focal->get_F();
  
  TorqueCalcSAM tocal(nmol);
  tocal.calc_all_tau(sys, focal);
  shared_ptr<vector<Pt> > to = tocal.get_tau();

  for (int i = 0; i < sys->get_n(); i++)
  {
    vector<Pt> fI = focal->get_all_fIk(i);
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      for (int d = 0; d < 3; d++)
      {
        if (k==0)
        {
          EXPECT_NEAR(for3[i][d]/(*fo)[i].get_cart(d), 1.0, preclim);
          EXPECT_NEAR(tor3[i][d]/(*to)[i].get_cart(d), 1.0, preclim);
        }
        if (fabs(forIk3[i][k][d]) > 1e-11)
          EXPECT_NEAR(forIk3[i][k][d]/fI[k].get_cart(d), 1.0, preclim);
      }
    }
//    cout << "This is force "<< setprecision(9) << fo[i].x()
//    << ", " << fo[i].y()<< ", " << fo[i].z()  << endl;
//    cout << "This is torque "<< setprecision(9) << to[i].x()
//    << ", " << to[i].y()<< ", " << to[i].z()  << endl;
  }
}


#endif /* EnergyForceUnitTest_h */

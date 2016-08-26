//
//  BDUnitTest.h
//  pbsam
//
//  Created by Lisa Felberg on 8/15/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef BDUnitTest_h
#define BDUnitTest_h

#include "BD.h"

/*
 Class for unit testing a BD step
 */
class BDStepUTest : public ::testing::Test
{
  protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}

  
};


TEST_F(BDStepUTest, three_mol_test)
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
  auto sys = make_shared<System>(mols);
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
  
  for (int i = 0; i < sys->get_n(); i++)
  {
//    vector<Pt> fI = focal->get_all_fIk(i);
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
//      for (int d = 0; d < 3; d++)
//      {
//        if (k==0)
//        {
//          EXPECT_NEAR(for3[i][d]/(*fo)[i].get_cart(d), 1.0, preclim);
//          EXPECT_NEAR(tor3[i][d]/(*to)[i].get_cart(d), 1.0, preclim);
//        }
//        if (fabs(forIk3[i][k][d]) > 1e-11)
//          EXPECT_NEAR(forIk3[i][k][d]/fI[k].get_cart(d), 1.0, preclim);
//      }
    }
    //    cout << "This is force "<< setprecision(9) << fo[i].x()
    //    << ", " << fo[i].y()<< ", " << fo[i].z()  << endl;
    //    cout << "This is torque "<< setprecision(9) << to[i].x()
    //    << ", " << to[i].y()<< ", " << to[i].z()  << endl;
  }
}

/*
 Class for unit testing a BD run
 */
class BDRunUTest : public ::testing::Test
{
  protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}
  
  
};


TEST_F(BDRunUTest, three_mol_test)
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
  auto sys = make_shared<System>(mols);
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
                       solvTest.get_IE(),
                       solvTest.get_interpol_list(), _expcons, pol);
  gsolvTest.solve(1e-16, 5);
  
  auto focal = make_shared<ForceCalc> (nmol, sys->get_all_Ik(),
                                       cst->get_dielectric_water(),
                                       SHCalcTest, BesselCal);
  focal->calc_all_f(solvTest.get_all_H(), solvTest.get_all_LHN(),
                    gsolvTest.get_gradH_all(),gsolvTest.get_gradLHN_all());
  shared_ptr<vector<Pt> > fo = focal->get_all_f();
  
  TorqueCalc tocal(nmol);
  tocal.calc_all_tau(sys, focal);
  shared_ptr<vector<Pt> > to = tocal.get_all_tau();
  
  for (int i = 0; i < sys->get_n(); i++)
  {
//    vector<Pt> fI = focal->get_all_fIk(i);
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      
//      for (int d = 0; d < 3; d++)
//      {
//        if (k==0)
//        {
//          EXPECT_NEAR(for3[i][d]/(*fo)[i].get_cart(d), 1.0, preclim);
//          EXPECT_NEAR(tor3[i][d]/(*to)[i].get_cart(d), 1.0, preclim);
//        }
//        if (fabs(forIk3[i][k][d]) > 1e-11)
//          EXPECT_NEAR(forIk3[i][k][d]/fI[k].get_cart(d), 1.0, preclim);
//      }
    }
    //    cout << "This is force "<< setprecision(9) << fo[i].x()
    //    << ", " << fo[i].y()<< ", " << fo[i].z()  << endl;
    //    cout << "This is torque "<< setprecision(9) << to[i].x()
    //    << ", " << to[i].y()<< ", " << to[i].z()  << endl;
  }
}

#endif /* BDUnitTest_h */

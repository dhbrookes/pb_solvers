//
//  EnergyForceUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 7/25/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef EnergyForceUnitTest_h
#define EnergyForceUnitTest_h

#include "PhysCalc.h"
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
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
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
    IEMatrix ieMatTest(0, sys->get_molecule(i),
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
  
  EnergyCalc ecal;
  vector<double> ene;
  ene = ecal.calc_all_energy(solvTest.get_all_H(), solvTest.get_all_LHN());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    EXPECT_NEAR(ene[i]/-0.00042138573, 1.0, preclim);
  }
}


TEST_F(EnergyUTest, three_mol_test)
{
  int pol(3), nmol(3);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
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
    IEMatrix ieMatTest(0, sys->get_molecule(i),
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
  
  EnergyCalc ecal;
  vector<double> ene;
  ene = ecal.calc_all_energy(solvTest.get_all_H(), solvTest.get_all_LHN());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
//    cout << "This is energy "<< setprecision(9) << ene[i] << endl;
    EXPECT_NEAR(ene[i]/en3[i], 1.0, preclim);
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
  
  vector<double> for3 = {};
  
};

TEST_F(ForceUTest, two_mol_test)
{
  int pol(3), nmol(2);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
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
    IEMatrix ieMatTest(0, sys->get_molecule(i),
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
//  solvTest.solve(1e-15, 200);
  
//  for (int i = 0; i < sys->get_n(); i++)
//  {
//    for (int k = 0; k < sys->get_Ns_i(i); k++)
//    {
//      solvTest.get_all_H()[i]->print_kmat(k);
//    }
//  }
  
  GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
                       solvTest.get_all_F(), solvTest.get_all_H(),
                       solvTest.get_IE(),
                       solvTest.get_interpol_list(), _expcons, pol);
  gsolvTest.solve(1e-16, 2);
  
  ForceCalc focal(SHCalcTest, BesselCal);
  vector<Pt> fo;
  fo = focal.calc_all_f(solvTest.get_all_H(), solvTest.get_all_LHN(),
                        gsolvTest.get_gradH_all(), gsolvTest.get_gradLHN_all());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    cout << "This is force "<< setprecision(9) << fo[i].x()
    << ", " << fo[i].y()<< ", " << fo[i].z()  << endl;
  }
}


TEST_F(ForceUTest, three_mol_test)
{
  int pol(3), nmol(3);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
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
    IEMatrix ieMatTest(0, sys->get_molecule(i),
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
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      solvTest.get_all_H()[i]->print_kmat(k);
    }
  }
  
//
//  GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
//                       solvTest.get_all_F(), solvTest.get_all_H(),
//                       solvTest.get_IE(),
//                       solvTest.get_interpol_list(), _expcons, pol);
//  gsolvTest.solve(1e-16, 72);
//  
//  ForceCalc focal(SHCalcTest, BesselCal);
//  vector<Pt> fo;
//  fo = focal.calc_all_f(solvTest.get_all_H(), solvTest.get_all_LHN(),
//                        gsolvTest.get_gradH_all(), gsolvTest.get_gradLHN_all());
//  
//  for (int i = 0; i < sys->get_n(); i++)
//  {
//    cout << "This is force "<< setprecision(9) << fo[i].x()
//    << ", " << fo[i].y()<< ", " << fo[i].z()  << endl;
////    EXPECT_NEAR(ene[i]/en3[i], 1.0, preclim);
//  }
}


#endif /* EnergyForceUnitTest_h */

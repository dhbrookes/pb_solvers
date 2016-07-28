//
//  ElectrostaticsUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 7/27/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef ElectrostaticsUnitTest_h
#define ElectrostaticsUnitTest_h

#include "Electrostatics.h"

/*
 Class for unit testing energy calculations
 */
class ElectroUTest : public ::testing::Test
{
  protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}
  
  vector<double> en3 = {-0.000883187183, -0.000457557568,-0.000425629615};
  
};

TEST_F(ElectroUTest, one_mol_test)
{
  int pol(3), nmol(1);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> (kT);
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
  string estart = test_dir_loc + "spol_test/test_0.00_p3.0.";
  vector<vector<string> > imat_loc(sys->get_n());
  string exp_str;
  
  vector<shared_ptr<HMatrix> > Hmat(sys->get_n());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    Hmat[i] = make_shared<HMatrix>(i, sys->get_Ns_i(i), pol, cst->get_kappa());
    imat_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
//      exp_loc[i][k] = estart+to_string(i)+"."+to_string(k)+ ".H.exp";
      exp_str = estart+to_string(k)+ ".H.exp";
      Hmat[i]->init_from_exp(exp_str, k);
    }
  }
  
  sys->write_to_pqr(test_dir_loc+"one_out.pqr");
  
  Electrostatic estat(Hmat, sys, SHCalcTest, BesselCal, cst, pol, 10);
  estat.print_3d_heat(test_dir_loc + "onemol_map.out");
  estat.print_dx(test_dir_loc + "one_mol.dx");
}

TEST_F(ElectroUTest, barnase_test)
{
  int pol(5), nmol(1);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<Molecule> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> (kT);
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
  string estart = test_dir_loc + "barnase_barstar/1BRS_chainA_0.05_p30.0.";
  vector<vector<string> > imat_loc(sys->get_n());
  string exp_str;
  
  vector<shared_ptr<HMatrix> > Hmat(sys->get_n());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    Hmat[i] = make_shared<HMatrix>(i, sys->get_Ns_i(i), pol, cst->get_kappa());
    imat_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
      //      exp_loc[i][k] = estart+to_string(i)+"."+to_string(k)+ ".H.exp";
      exp_str = estart+to_string(k)+ ".H.exp";
      Hmat[i]->init_from_exp(exp_str, k);
    }
  }
  
  sys->write_to_pqr(test_dir_loc+"barnase_out.pqr");
  
  Electrostatic estat(Hmat, sys, SHCalcTest, BesselCal, cst, pol, 50);
  estat.print_3d_heat(test_dir_loc + "barnase_map.out");
  estat.print_dx(test_dir_loc + "barnase_mol.dx");
}


TEST_F(ElectroUTest, three_mol_test)
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
  auto cst = make_shared<Constants> (kT);
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
  vector<vector<string > > exp_loc(sys->get_n());
  vector<shared_ptr<HMatrix> > Hmat(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    Hmat[i] = make_shared<HMatrix>(i, sys->get_Ns_i(i), pol, cst->get_kappa());
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
      exp_loc[i][k] = estart+to_string(i)+"."+to_string(k)+".H.exp";
      Hmat[i]->init_from_exp(exp_loc[i][k], k);
    }
  }
  
  sys->write_to_pqr(test_dir_loc+"three_out.pqr");
  Electrostatic estat(Hmat, sys, SHCalcTest, BesselCal, cst, pol, 10);
  estat.print_3d_heat(test_dir_loc + "threemol_map.out");
  estat.print_dx(test_dir_loc + "three_mol.dx");

}


#endif /* ElectrostaticsUnitTest_h */

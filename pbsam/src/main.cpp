/*
 main.cpp
 
 Main for PBSAM runs, electrostatics, dynamics and energyforce
 
 Copyright (c) 2015, Teresa Head-Gordon, Lisa Felberg, Enghui Yap, David Brookes
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of UC Berkeley nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "PBSAM.h"

using namespace std;

int main(int argc, const char * argv[])
{
//  string input_file = argv[1];
//  //  string input_file = "/Users/felb315/pb_solvers/pbsam/pbsam_test_files/energyforce_test/barnase_barstar/run.energyforce.hardrefs.inp";
//
//  PBSAM pbsam_run(input_file);
//  pbsam_run.run();
//  
//  return 0;
  string test_dir_loc = "/Users/davidbrookes/Projects/pb_solvers/pbsam/pbsam_test_files/gtest/";
//
//  int pol(3), nmol(2);
//  PQRFile pqr(test_dir_loc + "test_cged.pqr");
//  vector<shared_ptr<BaseMolecule> > mols;
//  for (int i=0; i<nmol; i++)
//    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
//                                            pqr.get_atom_pts(), pqr.get_radii(),
//                                            pqr.get_cg_centers(),
//                                            pqr.get_cg_radii()));
//  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
//  mols[1]->translate(Pt(10.71213542,-7.35779167,-0.15628125), 1e14);
//  auto sys = make_shared<System>(mols);
//  auto cst = make_shared<Constants> ();
//  cst->set_dielectric_water(80);
//  cst->set_dielectric_prot(4);
//  cst->set_salt_concentration(0.01);
//  cst->set_temp(298.15);
//  cst->set_kappa(0.0325628352);
//  
//  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
//  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
//  auto BesselCons = make_shared<BesselConstants> (2*pol);
//  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
//  auto _expcons = make_shared<ExpansionConstants> (pol);
//  
//  // Generate surface integrals
//  for (int i=0; i<nmol; i++)
//    IEMatrix ieMatTest(0, sys->get_moli(i),
//                       SHCalcTest, pol, _expcons, true, 0, true);
//  
//  string istart = test_dir_loc + "imat_test/imat.sp";
//  string estart = test_dir_loc + "spol_test/test_0.00_p3.0.";
//  vector<vector<string> > imat_loc(sys->get_n());
//  vector<vector<vector<string > > > exp_loc(sys->get_n());
//  
//  // Generate surface integrals
//  for (int i = 0; i < sys->get_n(); i++)
//  {
//    imat_loc[i].resize(sys->get_Ns_i(i));
//    exp_loc[i].resize(sys->get_Ns_i(i));
//    for (int k = 0; k < sys->get_Ns_i(i); k++)
//    {
//      exp_loc[i][k].resize(2);
//      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
//      exp_loc[i][k][0] = estart+to_string(k) + ".H.exp";
//      exp_loc[i][k][1] = estart+to_string(k) + ".F.exp";
//    }
//  }
//  
//  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
//                  true, true, imat_loc, exp_loc);
//  solvTest.solve(1e-15, 150);
  
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
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  //  cout << "This is cog of i " << sys->get_cogi(0).x() << ", " << sys->get_cogi(0).y() << ", " << sys->get_cogi(0).z() << endl;
  //  cout << "This is cog of i " << sys->get_cogi(1).x() << ", " << sys->get_cogi(1).y() << ", " << sys->get_cogi(1).z() << endl;
  //
  //  for (int j = 0; j < sys->get_n(); j++)
  //  {
  //    cout << "For MoleculeSAM  " << j << endl;
  //    for (int k = 0; k < sys->get_Ns_i(j); k++)
  //    {
  //      for (int d = 0; d < 3; d++)
  //        cout << sys->get_centerik(j,k).get_cart(d) << ",";
  //      cout << endl;
  //    }
  //  }
  
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
      imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
      exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+".H.exp";
      exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+".F.exp";
    }
  }
  
  Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                  true, true, imat_loc, exp_loc);
  solvTest.precalc_sh_lf_lh();
  solvTest.precalc_sh_numeric();
  GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
                       solvTest.get_all_F(), solvTest.get_all_H(),
                       solvTest.get_IE(), solvTest.get_interpol_list(),
                       solvTest.get_precalc_sh(), _expcons, pol);
  gsolvTest.solve(1e-16, 75);

}
//
//  SolverTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 6/22/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef SolverTest_h
#define SolverTest_h

#include "Solver.h"

class SolverTest
{
public:
  void RunSolverTest( )
  {
//    int pol(3), nmol(2);
//    PQRFile pqr(test_dir_loc + "test_cged.pqr");
//    vector<shared_ptr<Molecule> > mols;
//    for (int i=0; i<nmol; i++)
//      mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
//                                           pqr.get_atom_pts(), pqr.get_radii(),
//                                           pqr.get_cg_centers(),
//                                           pqr.get_cg_radii()));
//    
//    mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
//    mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
//    auto sys = make_shared<System>(mols);
//    auto cst = make_shared<Constants> ();
//    cst->set_dielectric_water(80);
//    cst->set_dielectric_prot(4);
//    cst->set_salt_concentration(0.01);
//    cst->set_temp(298.15);
//    cst->set_kappa(0.0325628352);
//    
//    auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
//    auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
//    auto BesselCons = make_shared<BesselConstants> (2*pol);
//    auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
//    auto _expcons = make_shared<ExpansionConstants> (pol);
//    
//    // Generate surface integrals
//    for (int i=0; i<nmol; i++)
//      IEMatrix ieMatTest(0, sys->get_molecule(i),
//                         SHCalcTest, pol, _expcons, true, 0, true);
//    
//    string istart = test_dir_loc + "imat_test/imat.sp";
//    string estart = test_dir_loc + "grad_test/mpol.";
//    vector<vector<string> > imat_loc(sys->get_n());
//    vector<vector<vector<string > > > exp_loc(sys->get_n());
//    
//    // Generate surface integrals
//    for (int i = 0; i < sys->get_n(); i++)
//    {
//      imat_loc[i].resize(sys->get_Ns_i(i));
//      exp_loc[i].resize(sys->get_Ns_i(i));
//      for (int k = 0; k < sys->get_Ns_i(i); k++)
//      {
//        exp_loc[i][k].resize(2);
//        imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
//        exp_loc[i][k][0] = estart+to_string(i)+"."+to_string(k)+".H.exp";
//        exp_loc[i][k][1] = estart+to_string(i)+"."+to_string(k)+".F.exp";
//      }
//    }
    
//    Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
//                    true, true, imat_loc, exp_loc);
//    GradSolver gsolvTest(sys, cst, SHCalcTest, BesselCal, solvTest.get_T(),
//                         solvTest.get_all_F(), solvTest.get_all_H(),
//                         solvTest.get_IE(), solvTest.get_interpol_list(),
//                         _expcons, pol);
//    gsolvTest.pre_compute_gradT_A();
//    gsolvTest.iter(0, 1);
  
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
  
};

#endif /* SolverTest_h */

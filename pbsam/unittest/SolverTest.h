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
    int pol(3), nmol(2);
    PQRFile pqr(test_dir_loc + "test_cged.pqr");
    vector<shared_ptr<Molecule> > mols;
    for (int i=0; i<nmol; i++)
      mols.push_back(make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                           pqr.get_atom_pts(), pqr.get_radii(),
                                           pqr.get_cg_centers(),
                                           pqr.get_cg_radii()));
    
    mols[0]->translate(Pt(-9.28786, -7.35779, -0.156281), 1e14);
    mols[1]->translate(Pt(10.71214, -7.35779, -0.156281), 1e14);
    auto sys = make_shared<System>(mols);
    
    cout << "This is cog of i " << sys->get_cogi(0).x() << ", " << sys->get_cogi(0).y() << ", " << sys->get_cogi(0).z() << endl;
    cout << "This is cog of i " << sys->get_cogi(1).x() << ", " << sys->get_cogi(1).y() << ", " << sys->get_cogi(1).z() << endl;
    auto cst = make_shared<Constants> ();
    cst->set_dielectric_water(80);
    cst->set_dielectric_prot(4);
    
    auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
    auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
    auto BesselCons = make_shared<BesselConstants> (2*pol);
    auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
    auto _expcons = make_shared<ExpansionConstants> (pol);
    
    string istart = test_dir_loc + "imat_test/imat.sp";
    string estart = test_dir_loc + "spol_test/test_0.00_p3.0.";
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
        exp_loc[i][k][0] = estart+to_string(k) + ".H.exp";
        exp_loc[i][k][1] = estart+to_string(k) + ".F.exp";
      }
    }
    
    Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol,
                    true, true, imat_loc, exp_loc);
  }
  
};

#endif /* SolverTest_h */

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
    
    int pol = 3;
    PQRFile pqr(test_dir_loc + "test_cged.pqr");
    auto myMol = make_shared<Molecule>(0, 0, "stat", pqr.get_charges(),
                                       pqr.get_atom_pts(), pqr.get_radii(),
                                       pqr.get_cg_centers(), pqr.get_cg_radii());
    auto cst = make_shared<Constants> ();
    cst->set_dielectric_water(80);
    cst->set_dielectric_prot(4);
    auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
    auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
    auto BesselCons = make_shared<BesselConstants> (2*pol);
    auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
    auto _expcons = make_shared<ExpansionConstants> (pol);
    
    // Generate surface integrals
    IEMatrix ieMatTest(0, myMol, SHCalcTest, pol, _expcons, true, 0, true);
    
    vector<shared_ptr<Molecule> > mols;
    mols.push_back(myMol);
    auto sys = make_shared<System>(mols);
    
    Solver solvTest( sys, cst, SHCalcTest, BesselCal, pol);
    solvTest.iter(0);
  }
  
};

#endif /* SolverTest_h */

//
//  SystemTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 06/01/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef SystemTest_h
#define SystemTest_h

#include "System.h"

using namespace std;

/*
 Class for testing spherical harmonics constants
 */
class SystemTest
{
public:
  void TestSystem( )
  {
  	 srand(1);
     PQRFile pqr(test_dir_loc + "test.pqr");
     MSMSFile surf_file (test_dir_loc + "test.vert");
     vector<Molecule> mols;
     mols.push_back(Molecule( 0, 0, "stat", pqr.get_charges(),
                       pqr.get_atom_pts(), pqr.get_radii(),
                       surf_file.get_sp(), surf_file.get_np(), 2.5));

  }

};



#endif


//
//  main.cpp
//  pb_solvers_test
//
//  Created by Marielle Soniat on 10/5/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#include <iostream>
#include <limits.h>
#include "gtest/gtest.h"

using namespace std;

double preclim = 1e-7;    //! precision limit
const int nvals = 10 ;      //! standard number of poles for testing
string test_dir_loc = "/Users/lfelberg/PBSAM/pb_solvers/pb_solvers_code/test/";

#include "ASolverUnitTest.h"
#include "BDUnitTest.h"
#include "BesselCalcUnitTest.h"
#include "ConstantsUnitTest.h"
#include "ElectrostaticsUnitTest.h"
#include "EnergyForceUnitTest.h"
#include "MyMatrixUnitTest.h"
#include "readutilUnitTest.h"
#include "ReExpCalcUnitTest.h"
#include "setupUnitTest.h"
#include "SHCalcUnitTest.h"
#include "SystemUnitTest.h"
#include "ThreeBodyUnitTest.h"
#include "utilUnitTest.h"

#include "tester.h"

int main(int argc, char * argv[])
{

  bool test = true;
  if (test)
  {
    cout << "Welcome to test suite" << endl;
    CTester allTest = CTester();
    allTest.unitTest( argc, argv );
  }
  
  vector<Molecule> mol;
  const int ml = 2;
  Pt pos[ml] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 6.0, 0.0)};
  for (int mi = 0; mi < ml; mi ++ )
  {
    int M = 3; vector<double> charges(M);
    vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 7.0; posCharges[0] = pos[mi];
    charges[1] = 7.0; posCharges[1] = pos[mi]+Pt(1,0,0);
    charges[2] = 7.0; posCharges[2] = pos[mi]+Pt(0,1,0);
    vdW[0]=0.0; vdW[1]=0.0; vdW[2]=0.0;
    Molecule molNew( "trans", 2.0, charges, posCharges, vdW, pos[mi], 0, 0.1);
    mol.push_back( molNew );
  }
  const int vals = 5;
  Constants const_;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver>( bCalcu, SHCalcu, sys, make_shared<Constants> (const_), vals);
  
  shared_ptr<TimeTerminate> term = make_shared<TimeTerminate>(10);
  BDRun BDTest( ASolvTest, term, false);
  BDTest.run();
//
  for (int step=0; step < 10; step ++)
  {
    //    cout <<setprecision(9)<< BDTest.get_system()->get_posij(2, 2).z() << ",";
    //    cout <<BDTest.get_system()->get_posij(0, 2).y() << ", " <<BDTest.get_system()->get_posij(0, 2).z() << " " << endl;
  }
  
  bool gtest = true;
  if (gtest)
  {
    cout << "Now for Google Tests." << endl ;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
  

  
}

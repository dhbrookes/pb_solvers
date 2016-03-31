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
string test_dir_loc = "/Users/lfelberg/PBSAM/pb_solvers/pbam/pbam_test_files/gtest/";

#include "ASolverUnitTest.h"
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
//#include "ThreeBodyUnitTest.h"
#include "utilUnitTest.h"

#include "BDUnitTest.h"

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
  
  bool gtest = true;
  if (gtest)
  {
    cout << "Now for Google Tests." << endl ;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
  

  
}

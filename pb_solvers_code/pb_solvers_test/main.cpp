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

double preclim = 1e-7;    //! precision limit
const int nvals = 10 ;      //! standard number of poles for testing

#include "ASolverUnitTest.h"
#include "BDUnitTest.h"
#include "BesselCalcUnitTest.h"
#include "ConstantsUnitTest.h"
#include "EnergyForceUnitTest.h"
#include "MyMatrixUnitTest.h"
#include "readutilUnitTest.h"
#include "ReExpCalcUnitTest.h"
#include "setupUnitTest.h"
#include "SHCalcUnitTest.h"
#include "SystemUnitTest.h"
#include "utilUnitTest.h"

#include "tester.h"

using namespace std;

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

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

#include "Constants.h"
#include "util.h"

const int nPol = 10 ;
double precLim = 1.0e-4;

#include "BesselCalcUnitTest.h"
#include "SHCalcUnitTest.h"

#include "tester.h"




using namespace std ;

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
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
}
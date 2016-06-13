/*
 main.cpp 
 
 Main for testing
 
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

#include <iostream>
#include <limits.h>
#include "gtest/gtest.h"

using namespace std;

double preclim = 1e-7;    //! precision limit
const int nvals = 10 ;      //! standard number of poles for testing
//string test_dir_loc = "../pbsam/pbsam_test_files/gtest/";
string test_dir_loc = "/Users/felb315/pb_solvers/pbsam/pbsam_test_files/gtest/";

#include "BesselCalcUnitTest.h"
#include "ConstantsUnitTest.h"
#include "MyMatrixUnitTest.h"
#include "readutilUnitTest.h"
#include "setupUnitTest.h"
#include "SHCalcUnitTest.h"
#include "SolverUnitTest.h"
#include "SystemUnitTest.h"
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
  
  bool gtest = true;
  if (gtest)
  {
    cout << "Now for Google Tests." << endl ;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
  

  
}

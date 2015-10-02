//
//  SHCalcTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/1/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef SHCalcTest_h
#define SHCalcTest_h

#include "SHCalc.h"

using namespace std;

/*
 Class for testing spherical harmonics constants
 */
class SHCalcConstantsTest
{
public:
  void TestSHCalcConstants( )
  {
    const int nPol = 10;
    double precLim = 1.0e-4;
    SHCalcConstants SHConstTest( nPol );
    assert( SHConstTest.get_n() == nPol ); // make sure numVals stores right
		
		// ensure proper calculation of (2*l-1)!!
    double doubleFactorial[10] = { 1.00000000e+00,   1.00000000e+00,
      3.00000000e+00, 1.50000000e+01,   1.05000000e+02,
      9.45000000e+02, 1.03950000e+04,   1.35135000e+05,
      2.02702500e+06, 3.44594250e+07 };
		
		// ensure propler calculation of shConsts
		double SHConst1[1]  = {};
		double SHConst5[1]  = {};
		double SHConst10[1] = {};
    
    for (int i = 0; i < nPol; i++) // check that our prefactors are right
    {
      assert( abs(doubleFactorial[i]-SHConstTest.get_dub_fac_val( i ))
                        < precLim);
      
    }
  }
}; // end SHCalcConstantsTest

#endif /* SHCalcTest_h */

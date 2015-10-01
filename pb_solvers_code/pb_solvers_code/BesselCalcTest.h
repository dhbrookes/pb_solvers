//
//  BesselCalcTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 9/30/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef BesselCalcTest_h
#define BesselCalcTest_h

#include "BesselCalc.h"

using namespace std;

/*
 Class for testing Bessel recursion constants
 */
class BesselConstantsTest
{
  
public:
  void TestBesselConstants( )
  {
    int nPol = 10;
		double preFactors[10] = {  };
    BesselConstants bConstTest = new BesselConstants( nPol );
    assert( bConstTest.get_n() == nPol );
    for (int i = 0; i < nPol; i++)
      assert( abs(preFactors[i]-bConstTest.get_const_val( i )) < 1.0e-4);    
  }  
};


/*
 Testing class for modified bessel functions (spherical and standard)
 */
class BesselCalcTest
{
protected:

  
public:
	
};


#endif /* BesselCalcTest_h */

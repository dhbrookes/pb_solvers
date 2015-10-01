//
//  BesselCalcTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 9/30/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef BesselCalcTest_h
#define BesselCalcTest_h

#include <stdio.h>
#include <vector>
#include <assert.h>

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
    double preFactors[10] = { 
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
  int                nPoles_;  // order of the Bessel function
  BesselConstants*  _consts_;  //constants used in recursion: Lotan 2006 eq3
  
public:
  
  /*
   Constructor if recursion constants have not already been calculated:
   */
  BesselCalc(int N, BesselConstants* _consts);
  
  /*
   Calculate the modified sphereical bessel functions I and K 
    (MBF of the first and second kind, respectively).
   Input is desired number of iterations an output is a vector 
    containing the calculated value at every iteration
   */
  const vector<double> calc_mbfI(const int n, const double z) const;
  const vector<double> calc_mbfK(const int n, const double z) const;
  
};


#endif /* BesselCalcTest_h */

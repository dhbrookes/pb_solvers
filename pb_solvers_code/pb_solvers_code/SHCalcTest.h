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
    const int vals = 10;
    double precLim = 1.0e-4;
    SHCalcConstants SHConstTest( vals );
    assert( SHConstTest.get_n() == vals );
    // make sure numVals stores right
		
		// ensure proper calculation of (2*l-1)!!
    double doubleFactorial[10] = { 1.00000000e+00,   1.00000000e+00,
      3.00000000e+00, 1.50000000e+01,   1.05000000e+02,
      9.45000000e+02, 1.03950000e+04,   1.35135000e+05,
      2.02702500e+06, 3.44594250e+07 };
		
		// ensure proper calculation of legendre Consts
    double LegConst1N5[6] = {1.8, 2.25,  3.0 ,  4.5 ,  9.0 , 0.0};
    double LegConst2N5[6] = {0.8, 1.25,  2.0 ,  3.5 ,  8.0 , 0.0};
    
    // ensure proper calculation of shConstant
    double SHConst[10] = { 1.00000000e+00, 9.53462589e-02, 9.17469804e-03, 8.99653161e-04,   9.08786929e-05,   9.57945535e-06, 1.07101567e-06,   1.29879727e-07,   1.76743922e-08, 2.86716502e-09 };

    
    // check that our prefactors are right
    // First check const1 and const2 of legendre polynomials
    assert( abs( 0.0 - SHConstTest.get_leg_consts1_val( 0, 0))
           < precLim);
    assert( abs( 0.0 - SHConstTest.get_leg_consts1_val( 0, 0))
           < precLim);
    
    for (int constIt = 0; constIt < 6; constIt++)
    {
      assert( abs( LegConst1N5[constIt] -
             SHConstTest.get_leg_consts1_val( 5, constIt)) < precLim);
      assert( abs( LegConst2N5[constIt] -
             SHConstTest.get_leg_consts2_val( 5, constIt)) < precLim);
    }
    
    // Check double factorial
    for (int i = 0; i < vals; i++)
    {
      assert( abs(doubleFactorial[i]-SHConstTest.get_dub_fac_val( i ))
             < precLim);
      assert( abs(SHConst[i]-SHConstTest.get_sh_consts_val( 10, i ))
             < precLim);
    }
    
  }
};


/*
 Class for testing spherical harmonics calculations,
 and legendre polynomials
 */
class SHCalcTest
{
public:
  void TestSHCalc( )
  {
    const int vals = 10;
    double precLim = 1.0e-4;
    SHCalcConstants SHConstTest( vals );
    SHCalc SHCalcTest( vals, &SHConstTest );
    
    // First test LEGENDRE
    // Testing \theta = 0, z = cos(\theta) = 1
    SHCalcTest.calc_sh( 0.0, 0.0 );
    assert( abs(SHCalcTest.get_legendre_result( 0, 0) - 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10, 0) - 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10,10) - 0.0)
           < precLim );
    
    // Testing \theta = \pi/3.0, z = cos(\theta) = 0.5
    double largeLeg = 1.55370279e+08;
    SHCalcTest.calc_sh( M_PI/3.0, 0.0 );
    assert( abs(SHCalcTest.get_legendre_result( 0, 0) - 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10, 0)
                + 1.88228607e-01) < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10,10)
                - largeLeg )/largeLeg < precLim );
    
    // Testing \theta = 2.0\pi/3.0, z = cos(\theta) = -0.5
    // Will be the same as z = 0.5
    SHCalcTest.calc_sh( 2.0*M_PI/3.0, 0.0 );
    assert( abs(SHCalcTest.get_legendre_result( 0, 0) - 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10, 0)
                + 1.88228607e-01) < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10,10)
                - largeLeg )/largeLeg < precLim );

    // Testing \theta = \pi, z = cos(\theta) = -1.0
    SHCalcTest.calc_sh( M_PI, 0.0 );
    assert( abs(SHCalcTest.get_legendre_result( 0, 0) - 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result( 5, 0) + 1.0)
           < precLim );
    assert( abs(SHCalcTest.get_legendre_result(10,10)
                - 0 ) < precLim );
    
    // Now test SPHERICAL HARMONICS
    // Testing \theta = 0, \phi = 0.0
    SHCalcTest.calc_sh( 0.0, 0.0 );
    assert( abs(SHCalcTest.get_result( 0, 0).real() - 1.0)
           < precLim ); // real component
    assert( abs(SHCalcTest.get_result( 0, 0).imag() - 0.0)
           < precLim ); // imag component
    assert( abs(SHCalcTest.get_result( 5, 5).real() - 0.0)
           < precLim ); // real component
    assert( abs(SHCalcTest.get_result( 5, 5).imag() - 0.0)
           < precLim ); // imag component
    
    // Testing \theta = 0.5, \phi = 0.5
    SHCalcTest.calc_sh( 0.5, 0.5 );
    assert( abs(SHCalcTest.get_result( 5, 0).real() + 0.16928726 )
           < precLim ); // real component
    assert( abs(SHCalcTest.get_result( 0, 0).imag() - 0.0)
           < precLim ); // imag component
    assert( abs(SHCalcTest.get_result( 5, 5).real() + 0.010066221)
           < precLim ); // real component
    assert( abs(SHCalcTest.get_result( 5, 5).imag() - 0.00751969)
           < precLim ); // imag component

  }
}; // end SHCalcTest

#endif /* SHCalcTest_h */

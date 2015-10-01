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
    const int nPol = 10;
		double precLim = 1.0e-4;
		BesselConstants bConstTest = BesselConstants( nPol );
		assert( bConstTest.get_n() == nPol ); // make sure numVals stores right
	
		double kPreFactors[10] = { -1.        ,
																0.33333333,
			                          0.06666667,
																0.02857143,
			                          0.01587302,
			                          0.01010101,
			                          0.00699301,
			                          0.00512821,
			                          0.00392157,
			                          0.00309598 };

    for (int i = 0; i < nPol; i++) // check that our prefactors are right
      assert( abs(kPreFactors[i]-bConstTest.get_kconst_val( i )) < precLim);
	}
};


/*
 Testing class for modified bessel functions (spherical and standard)
 */
class BesselCalcTest
{
public:
	void TestBesselCalc( )
	{
		const int nPol = 10;
		double precLim = 1.0e-4;
		BesselConstants * BConst = new BesselConstants( nPol );
		BesselCalc bConstTest = BesselCalc( nPol,  BConst );
		
		// for z = 1.0 and z = 10.0, calculated from python pbam_unit_test.py
		double i1[] = {1.17520119e+00,   1.10363832e+00,   1.07344305e+00,
		      1.05683451e+00,   1.04633846e+00,   1.03910889e+00,
		      1.03382786e+00,   1.02980185e+00,   1.02663135e+00,
			    1.02407006e+00};

		double i10[] = {1.10132329e+03,   2.97357289e+02,   1.20594900e+02,
		      6.18668362e+01,   3.69986800e+01,   2.46194746e+01,
		      1.77022637e+01,   1.34885612e+01,   1.07449415e+01,
		      8.86189182e+00};
		
		double k1[] = {1.00000000e+00,   2.00000000e+00,   2.33333333e+00,
		      2.46666667e+00,   2.53333333e+00,   2.57248677e+00,
		      2.59807600e+00,   2.61606542e+00,   2.62938888e+00,
			    2.63964796e+00 };

		double k10[] = {1.00000000e+00,   1.10000000e+01,   4.43333333e+01,
					1.17666667e+02,   2.44333333e+02,   4.31105820e+02,
					6.77907167e+02,   9.79379768e+02,   1.32702447e+03,
					1.71109497e+03};
		
		const vector<double> mBFI10 = bConstTest.calc_mbfI( nPol, 10.0 );
		const vector<double> mBFK10 = bConstTest.calc_mbfK( nPol, 10.0 );
		
		const vector<double> mBFI1 = bConstTest.calc_mbfI( nPol, 1.0 );
		const vector<double> mBFK1 = bConstTest.calc_mbfK( nPol, 1.0 );

		for (int besselIt = 0; besselIt < mBFI1.size(); besselIt++)
		{
			/*cout << " Program val " << mBFK1[besselIt] <<
			        " Python val " << k1[besselIt] << endl;*/
			assert( abs(mBFI1[besselIt] - i1[besselIt]) < precLim );
			assert( abs(mBFI10[besselIt] - i10[besselIt]) < precLim );
			assert( abs(mBFK1[besselIt] - k1[besselIt]) < precLim );
			assert( abs(mBFK10[besselIt] - k10[besselIt]) < precLim );

		}
	}
};


#endif /* BesselCalcTest_h */

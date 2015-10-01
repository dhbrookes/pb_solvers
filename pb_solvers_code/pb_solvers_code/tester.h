//
//  tester.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/1/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef tester_h
#define tester_h

#include <stdio.h>
#include <iostream>
#include "Constants.h"
#include "MyMatrix.h"
#include "ASolver.h"
#include "BesselCalc.h"
#include "BesselCalcTest.h"
#include "SHCalc.h"
#include "System.h"
#include "util.h"

/*
 Class for testing all Classes
 */
class CTester
{
public:
	
	int unitTest(int argc, const char * argv[])
	{
		BesselConstantsTest BCT = BesselConstantsTest();
		BCT.TestBesselConstants();
		
		BesselCalcTest BCalcT = BesselCalcTest();
		BCalcT.TestBesselCalc();
		
		cout << "Complete all tests" << endl;
		return 0;
	}

};

#endif /* tester_h */

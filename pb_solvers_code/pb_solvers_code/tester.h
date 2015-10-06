//
//  tester.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/1/15.
//  Copyright © 2015 LF. All rights reserved.
//

#ifndef tester_h
#define tester_h

#include <stdio.h>
#include <iostream>
#include "BesselCalcTest.h"
#include "ConstantsTest.h"
#include "SHCalcTest.h"

/*
 Class for testing all Classes
 */
class CTester
{
public:
	
	int unitTest(int argc, const char * argv[])
	{
		BesselConstantsTest BCT;
		BCT.TestBesselConstants();
		
		BesselCalcTest BCalcT;
		BCalcT.TestBesselCalc();
    
    ConstantsTest ConT;
    ConT.TestConstants();
		
		SHCalcConstantsTest SHConT;
		SHConT.TestSHCalcConstants();
    
    SHCalcTest SHCalcT;
    SHCalcT.TestSHCalc();
		
		cout << "Complete all tests" << endl;
		return 0;
	}

};

#endif /* tester_h */

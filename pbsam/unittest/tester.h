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
	
	int unitTest(int argc, char * argv[])
	{
    BesselConstantsTest BCT;
    BCT.TestBesselConstants();
    cout << "Complete Bessel constant test" << endl;
    
    BesselCalcTest BCalcT;
    BCalcT.TestBesselCalc();
    cout << "Complete Bessel calculation test" << endl;
    
//    ASolverTest ASolvT;
//    ASolvT.RunASolverTest();
//    cout << "Complete ASolv test" << endl;
    
    ConstantsTest ConT;
    ConT.TestConstants();
    cout << "Complete Constant test" << endl;
    
//    EnForTest EnForT;
//    EnForT.RunEnForTest();
//    cout << "Complete Energy and Force test" << endl;

//    ReExpTest ReR;
//    ReR.runReExTest();
//    cout << "Complete ReExpan test" << endl;
		
		SHCalcConstantsTest SHConT;
		SHConT.TestSHCalcConstants();
    cout << "Complete SH test" << endl;
    
    SHCalcTest SHCalcT;
    SHCalcT.TestSHCalc();
  
    cout << "Complete all tests" << endl;
    return 0;
  }

};

#endif /* tester_h */

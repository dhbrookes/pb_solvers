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
#include "ReExpCalcTest.h"
#include "SHCalcTest.h"
#include "SolverTest.h"
#include "SystemTest.h"

namespace pbsolvers
{

/*
 Class for testing all Classes
 */
class CTester
{
public:
  
  int unitTest(int argc, char * argv[])
  {
    SystemTest SysT;
    SysT.TestSystem();
    cout << "Complete system test" << endl;
    
    BesselConstantsTest BCT;
    BCT.TestBesselConstants();
    cout << "Complete Bessel constant test" << endl;
    
    BesselCalcTest BCalcT;
    BCalcT.TestBesselCalc();
    cout << "Complete Bessel calculation test" << endl;
    
    SolverTest SolvT;
    SolvT.RunSolverTest();
    cout << "Complete Solv test" << endl;
    
    ConstantsTest ConT;
    ConT.TestConstants();
    cout << "Complete Constant test" << endl;
    
    //    EnForTest EnForT;
    //    EnForT.RunEnForTest();
    //    cout << "Complete Energy and Force test" << endl;
    
    ReExpTest ReR;
    ReR.runReExTest();
    cout << "Complete ReExpan test" << endl;
    
    SHCalcConstantsTest SHConT;
    SHConT.TestSHCalcConstants();
    cout << "Complete SH test" << endl;
    
    SHCalcTest SHCalcT;
    SHCalcT.TestSHCalc();
    
    cout << "Complete all tests" << endl;
    return 0;
  }
  
};

} /* namespace pbsolvers */
#endif /* tester_h */

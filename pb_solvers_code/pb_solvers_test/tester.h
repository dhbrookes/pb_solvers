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
#include "ASolverTest.h"
#include "BesselCalcTest.h"
#include "EnergyForceTest.h"
#include "ConstantsTest.h"
#include "EnergyForceTest.h"
#include "ReExpCalcTest.h"
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
    
    ASolverTest ASolvT;
    ASolvT.RunASolverTest();
    cout << "Complete ASolv test" << endl;
    
    ConstantsTest ConT;
    ConT.TestConstants();
    cout << "Complete Constant test" << endl;
    
    EnForTest EnForT;
    EnForT.RunEnForTest();
    cout << "Complete Energy and Force test" << endl;

    ReExpTest ReR;
    ReR.runReExTest();
    cout << "Complete ReExpan test" << endl;
		
		SHCalcConstantsTest SHConT;
		SHConT.TestSHCalcConstants();
    cout << "Complete SH test" << endl;
    
    SHCalcTest SHCalcT;
    SHCalcT.TestSHCalc();
	
    vector<Molecule> mol_sing_;
    Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
      Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
      Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
    for (int molInd = 0; molInd < 9; molInd ++ )
    {
      int M = 3; vector<double> charges(M); vector<double> vdW(M);
      vector<Pt> posCharges(M);
      charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
      charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
      charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
      
      Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
      mol_sing_.push_back( molNew );
    }
    
    const int vals = 5;
    Constants const_( INTERNAL );
    shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
    shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
    shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
    shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
    shared_ptr<System> sys = make_shared<System>(mol_sing_);
    shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                          make_shared<Constants>
                                                          (const_), vals);
    ThreeBody threeBodTest( ASolvTest );
    threeBodTest.solveNmer(2);
  
    cout << "Complete all tests" << endl;
    return 0;
  }

};

#endif /* tester_h */

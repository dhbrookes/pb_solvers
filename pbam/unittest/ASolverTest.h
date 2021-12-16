//
//  ASolverTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/8/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ASolverTest_h
#define ASolverTest_h

#include "ASolver.h"

namespace pbsolvers
{

class ASolverTest
{
  
public :
  void RunASolverTest()
  {
    shared_ptr<Constants> const_ = make_shared<Constants>();
    vector< shared_ptr<BaseMolecule> > mol_;
    vector< shared_ptr<BaseMolecule> > mol_sing_;
    
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.5 ), Pt( 11.0, 6.9, 24.3 ) };
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, 0.0 ) };
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0] = cg[molInd]; posCharges[0] = cgPos[molInd]; vdW[0] = 0.0;
      
      shared_ptr<MoleculeAM> molNew = make_shared<MoleculeAM> ("stat",rd[molInd],charges,posCharges,vdW,pos[molInd],
                      molInd, 0);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0; posCharges[0] = cgPosSi[molInd];
      
      shared_ptr<MoleculeAM> molSing = make_shared<MoleculeAM> ( "stat", 10.0, charges, posCharges, vdW, molInd, 0);
      mol_sing_.push_back( molSing );
    }
    
    const int vals           = 10;
    shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
    shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
    shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
    shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
    shared_ptr<SystemAM> sys = make_shared<SystemAM>(mol_);
    
    ASolver ASolvTest (bCalcu, SHCalcu, sys, const_, vals);
    ASolvTest.solve_A( 1E-4 );
    ASolvTest.solve_gradA(1E-12);
  }
};


} /* namespace pbsolvers */

#endif /* ASolverTest_h */

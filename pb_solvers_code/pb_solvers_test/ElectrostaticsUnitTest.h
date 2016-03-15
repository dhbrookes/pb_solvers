//
//  ElectrostaticsUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ElectrostaticsUnitTest_h
#define ElectrostaticsUnitTest_h

#include "Electrostatics.h"

class ElecTest
{
  
  public :
  void RunElecTest()
  {
    Constants const_;
    vector< Molecule > mol_;

    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -8.0 ), Pt( 0,0,0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -7.0 ), Pt( 0,0,1) };
    double cg[2] = {-5.0, -5.0};
    double rd[2] = { 3.6,  3.6};

    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 3; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      charges[1]=cg[molInd]; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
      charges[2]=cg[molInd]; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);
      
      cout << setprecision(5) << posCharges[0].x() << " " << posCharges[0].y() << " " << posCharges[0].z() << " " << charges[0] << " " << 0.74 << endl;
      cout << setprecision(5) << posCharges[1].x() << " " << posCharges[1].y() << " " << posCharges[1].z() << " " << charges[1] << " " << 0.74 << endl;
      cout << setprecision(5) << posCharges[2].x() << " " << posCharges[2].y() << " " << posCharges[2].z() << " " << charges[2] << " " << 0.74 << endl;
      cout << setprecision(5) << pos[molInd].x() << " " << pos[molInd].y() << " " << pos[molInd].z() << " " << 0.0 << " " << rd[molInd] << endl;
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
    }
    
    const int vals = 5;
    BesselConstants bConsta( 2*vals );
    BesselCalc bCalcu( 2*vals, make_shared<BesselConstants>(bConsta) );
    SHCalcConstants SHConsta( 2*vals );
    SHCalc SHCalcu( 2*vals, make_shared<SHCalcConstants>(SHConsta) );
    System sys( mol_ );
    ReExpCoeffsConstants re_exp_consts (const_.get_kappa(),
                                        sys.get_lambda(), vals);
    
    shared_ptr<ASolver> ASolvTest = make_shared<ASolver>(make_shared<BesselCalc> (bCalcu),
                      make_shared<SHCalc> (SHCalcu),
                      make_shared<System> (sys),
                      make_shared<Constants> (const_), vals);
    ASolvTest->solve_A( 1E-12 ); ASolvTest->solve_gradA(1E-12);

    Electrostatic EstatTest( ASolvTest, 10);
    EstatTest.print_dx("/Users/lfelberg/Desktop/test.dx");
  }
};

#endif /* ElectrostaticsUnitTest_h */

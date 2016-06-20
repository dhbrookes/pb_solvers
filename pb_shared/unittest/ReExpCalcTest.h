//
//  ReExpCalcTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/13/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ReExpCalcTest_h
#define ReExpCalcTest_h

#include "ReExpCalc.h"

class ReExpTest 
{
  public :
  
  ReExpTest() : vals_(nvals)   {  }
  
  void runReExTest()
  {
    Constants Cst;
    Pt testPt = Pt( 0.0, 0.0, 5.0);
//    SHCalcConstants shCon( 2*nvals );
    SHCalc shCalc( 2*nvals, make_shared<SHCalcConstants>(2*nvals) );
    shCalc.calc_sh( testPt.theta(), testPt.phi());
    
//    BesselConstants bCon( 2*nvals );
    BesselCalc      bCal( 2*nvals, make_shared<BesselConstants>(2*nvals) );
    
    MyMatrix<cmplx> shMat = shCalc.get_full_result();
    double kap            = Cst.get_kappa();
    double lambda         = 5.0;
    
    vector<double> besselK = bCal.calc_mbfK(2*nvals, kap*testPt.r());
    
    ReExpCoeffs ReExpTest( nvals, testPt, shMat, besselK,
                          make_shared<ReExpCoeffsConstants> (kap,lambda,nvals),
                          kap, {lambda}, true );
    
  }
  
  int get_vals()                { return vals_; }
  
  protected :
  
  int vals_;
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};



#endif /* ReExpCalcTest_h */

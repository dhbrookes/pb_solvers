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
    SHCalcConstants shCon( 2*nvals );
    SHCalc shCalc( 2*nvals, shCon );
    shCalc.calc_sh( testPt.theta(), testPt.phi());
    
    BesselConstants bCon( 2*nvals );
    BesselCalc      bCal( 2*nvals, bCon );
    
    MyMatrix<cmplx> shMat = shCalc.get_full_result();
    double kap            = Cst.get_kappa();
    double lambda         = 5.0;
    ReExpCoeffsConstants reExpConsts (kap, lambda, nvals);
    
    ReExpCoeffs ReExpTest( testPt, nvals, shMat, bCal, reExpConsts,
                          kap, lambda );
    
  }
  
  int get_vals()                { return vals_; }
  
  protected :
  
  int vals_;
  
  virtual void SetUp()     { }
  virtual void TearDown()  { }
};



#endif /* ReExpCalcTest_h */

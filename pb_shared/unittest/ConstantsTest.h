//
//  ConstantsTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/6/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef ConstantsTest_h
#define ConstantsTest_h

#include "Constants.h"
using namespace std;

namespace pbsolvers
{

/*
 Class for testing system constants class
 */
class ConstantsTest
{
public:
  void TestConstants( )
  {
    double precLim = 1.0e-4;
    Constants ConstTest;
    
    // With system defaults, testing kappa calc
    assert( abs(ConstTest.get_kappa() - 0.03030729144) < precLim );
    
    // check Kappa after changing \( eps_s \)
    ConstTest.set_dielectric_water( 40.0 );
    assert( abs(ConstTest.get_kappa() - 0.04232182927) < precLim );
    
    // check Kappa after changing Salt concentration \( Molar \)
    ConstTest.set_salt_concentration( 0.05 );
    assert( abs(ConstTest.get_kappa() - 0.09463448717) < precLim );
    
    // check Kappa after changing temperature [ Kelvin ]
    ConstTest.set_temp( 298.0 );
    assert( abs(ConstTest.get_kappa() - 0.1029979673) < precLim );
    assert( abs(ConstTest.get_kbt()   - 4.11436084E-21) < precLim );
    double iKt = 2.430511175E20; // large number
    assert( abs(ConstTest.get_ikbt()  - iKt)/iKt < precLim );
    
    // check unit conversion
    assert( abs(ConstTest.convert_int_to_kcal_mol( 1.0 ) - 332.061203)
           < precLim );
    assert( abs(ConstTest.convert_int_to_jmol( 1.0 ) - 1389344.0722)
           < precLim );
    assert( abs(ConstTest.convert_int_to_kT( 1.0 ) - 560.73826468)
           < precLim );
  }
  
};

} /* namespace pbsolvers */
#endif /* ConstantsTest_h */

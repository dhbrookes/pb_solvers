//
//  ConstantsUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/6/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef ConstantsUnitTest_h
#define ConstantsUnitTest_h

#include "Constants.h"
double preclim = 1.0e-4;

Constants ConstUTest;

class ConstantsUTest : public ::testing::Test
{
  protected :
  
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(ConstantsUTest, settingAndCalc)
{
  EXPECT_NEAR( ConstUTest.get_kappa(), 0.03030729144, preclim);
  
  // check Kappa after changing \( eps_s \)
  ConstUTest.set_dielectric_water( 40.0 );
  EXPECT_NEAR( ConstUTest.get_kappa(), 0.04232182927, preclim);
  
  // check Kappa after changing Salt concentration \( Molar \)
  ConstUTest.set_salt_concentration( 0.05 );
  EXPECT_NEAR( ConstUTest.get_kappa(), 0.09463448717, preclim);
  
  // check Kappa after changing temperature [ Kelvin ]
  ConstUTest.set_temp( 298.0 );
  EXPECT_NEAR( ConstUTest.get_kappa(), 0.1029979673, preclim);
  EXPECT_NEAR( ConstUTest.get_kbt(), 4.11436084E-21, preclim);
  double iKt = 2.430511175E20; // large number
  EXPECT_NEAR( ConstUTest.get_ikbt()/iKt, 1.0, preclim);
}

TEST_F(ConstantsUTest, unitConv)
{
  // check unit conversion
  EXPECT_NEAR( ConstUTest.convert_int_to_kcal_mol( 1.0 ), 332.061203, preclim);
  EXPECT_NEAR( ConstUTest.convert_int_to_jmol( 1.0 ), 1389344.0722, preclim);
  EXPECT_NEAR( ConstUTest.convert_int_to_kT( 1.0 ), 560.73826468, preclim);
  
}

#endif /* ConstantsUnitTest_h */

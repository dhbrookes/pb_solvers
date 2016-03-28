//
//  setupUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef setupUnitTest_h
#define setupUnitTest_h

#include "setup.h"

/*
 Class to test setting up details for the PBAM system
 */
class SetupUTest : public ::testing::Test
{
public :
  
protected :

};

TEST_F(SetupUTest, checkOpen)
{
  string path = test_dir_loc + "none.inp";
  try
  {
    class Setup setTest(path);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(SetupUTest, readInp)
{
  string inp = test_dir_loc + "run.1.inp";
  class Setup setTest(inp);
 
  ASSERT_EQ("dynamics", setTest.getRunType());
  ASSERT_EQ("test1", setTest.getRunName());

  ASSERT_EQ( 353, setTest.getTemp());
  ASSERT_EQ( 0.0, setTest.getSaltConc());
  
  ASSERT_EQ( 1, setTest.getIDiel());
  ASSERT_EQ( 8, setTest.getSDiel());
  
  ASSERT_EQ( 3, setTest.getNType());
  
  ASSERT_EQ(  0.000, setTest.getDtr(0));
  ASSERT_EQ(  0.000, setTest.getDrot(0));
  ASSERT_EQ(  0.000, setTest.getDtr(1));
  ASSERT_EQ( 9.2e-5, setTest.getDrot(1));
  ASSERT_EQ(  0.811, setTest.getDtr(2));
  ASSERT_EQ( 0.5e-5, setTest.getDrot(2));
  
  ASSERT_EQ( "stat", setTest.getTypeNDef(0));
  ASSERT_EQ(  "rot", setTest.getTypeNDef(1));
  ASSERT_EQ( "move", setTest.getTypeNDef(2));
  
  ASSERT_EQ( "barnase.pqr", setTest.getTypeNPQR(0));
  ASSERT_EQ( "2charge.pqr", setTest.getTypeNPQR(1));
  ASSERT_EQ( "2charge.pqr", setTest.getTypeNPQR(2));
  
  ASSERT_EQ(  "lisa.xyz", setTest.getTypeNXYZ(0));
  ASSERT_EQ( "other.xyz", setTest.getTypeNXYZ(1));
  ASSERT_EQ(   "trd.xyz", setTest.getTypeNXYZ(2));
  
}


#endif /* setupUnitTest_h */

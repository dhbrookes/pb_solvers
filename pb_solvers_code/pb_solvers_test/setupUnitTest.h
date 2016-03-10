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
  // This will change according to where you're running tests!
  string test_dir_ = "/Users/lfelberg/PBSAM/pb_solvers/pb_solvers_code/test/";
  
};

TEST_F(SetupUTest, checkOpen)
{
  string path = test_dir_ + "none.inp";
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
  string inp = test_dir_ + "run.1.inp";
  class Setup setTest(inp);
 
//  ASSERT_EQ(true, PQRtest.get_cg());

}


#endif /* setupUnitTest_h */

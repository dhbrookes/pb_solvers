//
//  SystemUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/9/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef SystemUnitTest_h
#define SystemUnitTest_h

#include "System.h"

class MoleculeUTest : public ::testing::Test
{
  public :
  
  protected :
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};


TEST_F(MoleculeUTest, checkVals)
{
  
}


class SystemUTest : public ::testing::Test
{
  public :
  
  protected :
  
  virtual void SetUp() {}
  virtual void TearDown() {}
  
};


TEST_F(SystemUTest, checkVals)
{
  // EXPECT_NEAR( testMat3(1,1), 140.0, preclim);
}


#endif /* SystemUnitTest_h */

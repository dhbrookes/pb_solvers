//
//  main.cpp
//  Tests
//
//  Created by Marielle Soniat on 10/5/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#include <iostream>

#include <limits.h>
#include "gtest/gtest.h"

#include "Constants.h"
#include "util.h"

#include "BesselCalcUnitTest.h"
#include "SHCalcUnitTest.h"

using namespace std ;

int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

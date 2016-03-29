#include <limits.h>
#include "sample1.h"
#include "gtest/gtest.h"


TEST(add_ints, sample)
{
  EXPECT_EQ(7, add_ints(3,4)) << "This should pass." ;
  EXPECT_EQ(1, add_ints(0,1));
  EXPECT_EQ(7, add_ints(4,4)) << "Example of gtest failure." ;
}

#include "BesselCalc.h"

#include <limits.h>
#include "gtest/gtest.h"

TEST(BesselConstants, first10)
{
  const int nPol = 10 ;
  double precLim = 1.0e-4;
  BesselConstants bConstTest = BesselConstants( nPol );
  ASSERT_EQ( nPol , bConstTest.get_n() ) ;
}


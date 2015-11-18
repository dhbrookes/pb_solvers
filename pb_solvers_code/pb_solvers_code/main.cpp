//
//  main.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/24/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include <iostream>
#include "Constants.h"
#include "MyMatrix.h"
#include "ASolver.h"
#include "BesselCalc.h"
#include "SHCalc.h"
#include "System.h"
#include "util.h"
#include <memory>
#include "readutil.h"

using namespace std;

int main(int argc, const char * argv[])
{
  vector<double> av (5, 1);
  MyVector<double> a (av);
  cout << a[2] << endl;
}

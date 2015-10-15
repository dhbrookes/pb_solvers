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


int main(int argc, const char * argv[])
{
	
  cout << "Hello I build and run" << endl;
  
  vector<MyMatrix<double> > a (2, (2, 2));
  MyVector<MyMatrix<double> >* my_a = new MyVector<MyMatrix<double> >(a);
  MyMatrix<double>* b = &my_a->operator[](0);
  b->set_val(0, 0, 5);
  
	return 0;
}

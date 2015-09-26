//
//  BesselCalc.hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef BesselCalc_h
#define BesselCalc_h

#include <stdio.h>
#include <vector>

using namespace std;

/*
 Calculator class for modified bessel functions (spherical and standard)
 */
class BesselCalc
{
protected:
    
    int N_;  // order of the Bessel function
    vector<double> recConsts_;  //constants used in recursion: Lotan 2006 eq3
    
public:
    
    /*
     Constructor initializes constants
     */
    BesselCalc(int N);
    
    /*
     Calculate the modified sphereical bessel functions I and K (MBF of the first and second kind, respectively).
     Input is desired number of iterations an output is a vector containing the calculated value at every iteration
     */
    vector<double> calc_mbfI(int num_iter, double z);
    vector<double> calc_mbfK(int num_iter, double z);
    
};

#endif /* BesselCalc_hpp */

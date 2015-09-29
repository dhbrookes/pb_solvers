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
#include <assert.h>

using namespace std;


/*
 Class for storing constants that can be used for multiple runs of bessel calcualtions
 */
class BesselConstants
{
protected:
    int nPoles_;
    vector<double> recConsts_;  // recursion constants
    
public:
    
    BesselConstants(const int N);
    
    const int get_n() const             { return nPoles_; }
    const double get_const_val(int i)   { return recConsts_[i]; }
    
};


/*
 Calculator class for modified bessel functions (spherical and standard)
 */
class BesselCalc
{
protected:
    
    int                 nPoles_;  // order of the Bessel function
    BesselConstants*    _consts_;  //constants used in recursion: Lotan 2006 eq3
    
public:
    
    /*
     Constructuro if recursion constants have not already been calculated:
     */
    BesselCalc(int N, BesselConstants* _consts);
    
    /*
     Calculate the modified sphereical bessel functions I and K (MBF of the first and second kind, respectively).
     Input is desired number of iterations an output is a vector containing the calculated value at every iteration
     */
    const vector<double> calc_mbfI(const int n, const double z) const;
    const vector<double> calc_mbfK(const int n, const double z) const;
    
};


#endif /* BesselCalc_hpp */

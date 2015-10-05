//
//  BesselCalc.h
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef BesselCalc_h
#define BesselCalc_h

#include <stdio.h>
#include "Constants.h"
#include <vector>
#include <assert.h>

using namespace std;

/*
 Class for storing constants that can be used for multiple runs of 
  bessel calculations
 */
class BesselConstants
{
protected:
  int numVals_;
  vector<double> kConsts_;  // recursion constants for k
  
public:
    
    BesselConstants(const int N=Constants::MAX_NUM_POLES);
    
    const int get_n() const                 { return numVals_; }
    const double get_kconst_val(int i)      { return kConsts_[i]; }
    const double get_iconst_val(int i)      { return iConsts_[i]; }
    
};


/*
 Calculator class for modified bessel functions (spherical and standard)
 */
class BesselCalc
{
protected:
    
    int                 numVals_;  // order of the Bessel function
    BesselConstants*    _consts_;  // constants used in recursion: Lotan 2006 eq3
    
public:
    
    BesselCalc();
  
    BesselCalc(int N, BesselConstants* _consts);
    
    BesselCalc(const BesselCalc& other);
    
    virtual ~BesselCalc();
    
    BesselCalc& operator=(const BesselCalc& other);
  
  /*
   Calculate the modified sphereical bessel functions I and K 
  (MBF of the first and second kind, respectively).
   Input is desired number of iterations an output is a vector 
  containing the calculated value at every iteration
   */
  const vector<double> calc_mbfI(const int n, const double z) const;
  const vector<double> calc_mbfK(const int n, const double z) const;
  
};



#endif /* BesselCalc_h */

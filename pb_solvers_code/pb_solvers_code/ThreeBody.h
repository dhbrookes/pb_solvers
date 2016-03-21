//
//  ThreeBody.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/9/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ThreeBody_h
#define ThreeBody_h

#include "EnergyForce.h"

/*
 This class is designed to compute the vector A defined in Equation 22
 Lotan 2006, page 544
 */
class ThreeBody
{
protected:
  int N_; // number of molecules in the system
  int p_; // number of poles
  
  double cutoffTBD_;  // distance for cutoff of tbd approx
  
  vector<vector<int > > dimer_;   // list of all pairs by their index #
  vector<vector<int > > trimer_;  // list of all triplets by their index #
  
  vector<vector< double > > energy_di_;
  vector<vector< double > > energy_tri_;
  
  vector<vector< Pt > > force_di_;
  vector<vector< Pt > > force_tri_;
  
  shared_ptr<BesselCalc>      _besselCalc_;
  shared_ptr<System>          _sys_;  // system data (radii, charges, etc.)
  shared_ptr<SHCalc>          _shCalc_;
  shared_ptr<Constants>       _consts_;
  
private:
  ThreeBody();
  
  void generatePairsTrips();
  
  // Solve the N body problem, only 2 or 3 right now
  void solveNmer( int num);
  
  shared_ptr<System> make_subsystem(vector<int> mol_idx);
  
};

#endif /* ThreeBody_h */

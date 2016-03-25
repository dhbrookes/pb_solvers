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
  
  Units unit_; // String value of units
  
  vector<vector<int > > dimer_;   // list of all pairs by their index #
  vector<vector<int > > trimer_;  // list of all triplets by their index #
  
  vector<vector< double > > energy_di_;
  vector<vector< double > > energy_tri_;
  
  vector<vector< Pt > > force_di_;
  vector<vector< Pt > > force_tri_;
  
  vector<vector< Pt > > torque_di_;
  vector<vector< Pt > > torque_tri_;
  
  vector<double> energy_approx_;
  vector< Pt >   force_approx_;
  vector< Pt >   torque_approx_;
  
  shared_ptr<BesselCalc>      _besselCalc_;
  shared_ptr<System>          _sys_;  // system data (radii, charges, etc.)
  shared_ptr<SHCalc>          _shCalc_;
  shared_ptr<Constants>       _consts_;
  
  shared_ptr<System> make_subsystem(vector<int> mol_idx);
  
  int find_di( int i, int j);
  void generatePairsTrips();
  
public:
  ThreeBody(shared_ptr<ASolver> _asolver, Units unt = INTERNAL,
            double cutoff = 1e48 );
  
  // Solve the N body problem, only 2 or 3 right now
  void solveNmer( int num, double preclim = 1e-4);
  
  void printNmer( int num, string outfile);
  
  void calcTBDEnForTor( );
  void calcTwoBDEnForTor( );
  
  void printTBDEnForTor( vector<string> outfile );
  
  vector<vector<int > > getDimers()  { return dimer_; }
  vector<vector<int > > getTrimers() { return trimer_; }
  
  vector<vector<double > > getDiEn() { return energy_di_; }
  vector<vector<double > > getTrEn() { return energy_tri_; }
  
  vector<vector<Pt > > getDiFo()     { return force_di_; }
  vector<vector<Pt > > getTrFo()     { return force_tri_; }
  
  vector<vector<Pt > > getDiTo()     { return torque_di_; }
  vector<vector<Pt > > getTrTo()     { return torque_tri_; }
  
  vector<double> get_energy_approx() { return energy_approx_; }
  vector<Pt> get_force_approx()      { return force_approx_; }
  vector<Pt> get_torque_approx()     { return torque_approx_; }
  
  double get_energyi_approx( int i) { return energy_approx_[i]; }
  Pt get_forcei_approx( int i)      { return force_approx_[i]; }
  Pt get_torquei_approx( int i)     { return torque_approx_[i]; }
  
};

#endif /* ThreeBody_h */

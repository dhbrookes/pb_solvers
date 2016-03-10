//
//  Electrostatics.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef Electrostatics_h
#define Electrostatics_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "ASolver.h"

/*
 Class for printing out electrostatics of system
 */
class Electrostatic
{
protected:
  int p_; // Npoles
  
  vector<double> range_min_;  // Origin of grid in each dim
  vector<double> range_max_;  // Origin of grid in each dim
  vector<int> npts_;   // number of grid pts in each dimension
  vector<double> step_;  // step of grid in each dimension

  vector < vector < vector<double > > > esp_; // vector of ESP values
  
  VecOfMats<cmplx>::type A_;
  shared_ptr<System> _sys_;
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr<BesselCalc> _bCalc_;
  
  void find_range();
  void find_bins();
  
  void compute_pot();
  double compute_pot_at( double x, double y, double z);
  
  MyMatrix<cmplx> get_local_exp( Pt dist );
  
  double lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V, int p);
  
public:
  Electrostatic(VecOfMats<cmplx>::type A, System sys,
                SHCalc shCalc, BesselCalc bCalc, int p, int npts = 150);
  
  // print APBS file
  void print_dx(string ifname);
  
  void print_grid(string dim, double value);
  
};


#endif /* Electrostatics_h */

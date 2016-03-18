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

enum Axis {Xdim, Ydim, Zdim};

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

  vector<vector<vector<double > > > esp_; // vector of ESP values
  
  shared_ptr<VecOfMats<cmplx>::type> _A_;
  shared_ptr<System> _sys_;
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr<BesselCalc> _bCalc_;
  shared_ptr<Constants> _consts_;
  
  void find_range();
  void find_bins();
  
  void compute_pot();
  double compute_pot_at( double x, double y, double z);
  
  MyMatrix<cmplx> get_local_exp( Pt dist );
  
  double lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V, int p);
  
public:
  Electrostatic(shared_ptr<VecOfMats<cmplx>::type> _A, shared_ptr<System> _sys,
                shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
                shared_ptr<Constants> _consts,
                int p, int npts = 150);
  
  Electrostatic(shared_ptr<ASolver> _asolv, int npts=150);
  
  // print APBS file
  void print_dx(string ifname);
  
  void print_grid(Axis axis, double value, string fname);
  
};


#endif /* Electrostatics_h */

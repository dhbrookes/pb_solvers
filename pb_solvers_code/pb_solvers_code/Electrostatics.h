//
//  Electrostatics.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef Electrostatics_h
#define Electrostatics_h

#include "ASolver.h"

enum Axis {Xdim, Ydim, Zdim};

/*
 Class for printing out electrostatics of system
 */
class Electrostatic
{
protected:
  int p_; // Npoles
  double units_; // A conversion factor to user desired units
  
  double pot_min_; // A minimum value of the pot
  double pot_max_; // A max value of the potential
  
  vector<double> range_min_;  // Origin of grid in each dim
  vector<double> range_max_;  // Origin of grid in each dim
  vector<int> npts_;   // number of grid pts in each dimension
  vector<double> step_;  // step of grid in each dimension

  vector<vector<vector<double > > > esp_; // vector of ESP values
  vector<vector<double > > grid_;  // 2D cross section of ESP
  
  shared_ptr<VecOfMats<cmplx>::type> _A_;
  shared_ptr<System> _sys_;
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr<BesselCalc> _bCalc_;
  shared_ptr<Constants> _consts_;
  
  void find_range();
  void find_bins();
  
  void compute_units();
  
  void compute_pot();
  double compute_pot_at( Pt point );
  
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
  
  // print Grid file, given an axis and a value on that axis
  void print_grid(string axis, double value, string fname);
  
  // return potential grid
  vector<vector<vector<double > > > get_potential() { return esp_; }
  vector<vector<double > > get_pot2d()              { return grid_; }
  
  vector<double> get_mins()  { return range_min_; }
  vector<double> get_maxs()  { return range_max_; }
  vector<int> get_npts()     { return npts_; }
  vector<double> get_bins()  { return step_; }
  
};


#endif /* Electrostatics_h */

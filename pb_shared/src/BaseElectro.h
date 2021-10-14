//
//  BaseElectro.h
//  pbsam_xcode
//
//  Created by David Brookes on 8/31/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BaseElectro_h
#define BaseElectro_h

#include <stdio.h>

//#include "BasePhysCalc.h"
#include <time.h>

#ifdef __OMP
#include <omp.h>
#endif

#include <vector>
#include <string>
#include "BesselCalc.h"
#include "SHCalc.h"
#include "BaseSys.h"

using namespace std;

// Exception class to ensure that
class ValueOutOfRange: public exception
{
protected:
  string ax_;
  double value_;
  double range_;
  mutable string ss_;
  
public:
  ValueOutOfRange( string ax, double val, double ran)
  :ax_(ax), value_(val), range_(ran), ss_()
  {
  }
  
  virtual const char* what() const throw()
  {
    ss_ = ax_ + " value " + to_string(value_)+ " out of range. It is";
    if (value_ < range_)
    ss_ += " less than ";
    else
    ss_ += " greater than ";
    ss_ += to_string(range_);
    return ss_.c_str();
  }
};

/*
 Class for printing out electrostatics of system
 */
class BaseElectro
{
protected:
  int p_; // Npoles
  double units_; // A conversion factor to user desired units
  
  double pot_min_; // A minimum value of the pot
  double pot_max_; // A max value of the potential
  
  double lam_; // Average radius of MoleculeSAMs in system
  
  vector<double> range_min_;  // Origin of grid in each dim
  vector<double> range_max_;  // Origin of grid in each dim
  vector<int> npts_;   // number of grid pts in each dimension
  vector<double> step_;  // step of grid in each dimension
  
  vector<vector<vector<double > > > esp_; // vector of ESP values
  vector<vector<double > > grid_;  // 2D cross section of ESP
  
//  vector<shared_ptr<HMatrix> > _H_;
  shared_ptr<BaseSystem> _sys_;
  shared_ptr<SHCalc> _shCalc_;
  shared_ptr<BesselCalc> _bCalc_;
  shared_ptr<Constants> _consts_;
  
  void find_range();
  void find_bins();
  
  void compute_units();
  
  void compute_pot();
  
  // main method to be overridden in sub classes
  virtual double compute_pot_at( Pt point )=0;
  
  MyMatrix<cmplx> get_local_exp( Pt dist, double lambda );
  
  double lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V, int p);
  
public:
  BaseElectro(shared_ptr<BaseSystem> _sys,
              shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
              shared_ptr<Constants> _consts,
              int p, int npts = 150);
  
  
  // print APBS file
  void print_dx(string ifname);
  
  // print out 3D heatmap data for surface of each sphere
  void print_3d_heat( string td_name );
  
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


#endif /* BaseElectro_h*/

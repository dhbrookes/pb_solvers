//
//  BD.h
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BD_h
#define BD_h

#include <stdio.h>
#include "System.h"
#include "util.h"
#include <random>
#include <memory>
#include "ASolver.h"
#include "EnergyForce.h"

/*
 Base class for implementing termination conditions in BD
 */
class BaseTerminate
{
public:
  BaseTerminate() { }
  
  virtual const bool is_terminated(shared_ptr<System> _sys) const
  {
    return false;
  }
};

/*
 Class for time based termination
 */
class TimeTerminate : public BaseTerminate
{
protected:
  double endTime_; //termination time
  
public:
  TimeTerminate(double end_time)
  :BaseTerminate(), endTime_(end_time)
  {
  }
  
  const bool is_terminated(shared_ptr<System> _sys) const
  {
    bool done = false;
    if (_sys->get_time() >= endTime_) done = true;
    return done;
  }
};

enum CoordType { X, Y, Z, R };
enum BoundaryType { EQ, LEQ, GEQ };

/*
 Class for coordinate based termination. This terminates based on whether
 the specified molecule satisfies the BoundaryType condition on the CoordType
 with the given boundary value.
 */
class CoordTerminate : public BaseTerminate
{
protected:
  double boundaryVal_;
  int molIdx_;
  CoordType coordType_;
  BoundaryType boundType_;
  
public:
  CoordTerminate(int mol_idx, CoordType coord_type,
                 BoundaryType bound_type, double bound_val)
  :BaseTerminate(), molIdx_(mol_idx), coordType_(coord_type),
  boundType_(bound_type), boundaryVal_(bound_val)
  {
  }
  
  const bool is_terminated(shared_ptr<System> _sys) const
  {
    bool done;
    Pt mol_coord = _sys->get_centeri(molIdx_);
    double test_val;
    if (coordType_ == X)      test_val = mol_coord.x();
    else if (coordType_ == Y) test_val = mol_coord.y();
    else if (coordType_ == Z) test_val = mol_coord.z();
    else                      test_val = mol_coord.norm();
    
    if (boundType_ == EQ)
      test_val == boundaryVal_ ? done = true: done = false;
    else if (boundType_ == LEQ)
      test_val < boundaryVal_ ? done = true: done = false;
    else if (boundType_ == GEQ)
      test_val > boundaryVal_ ? done = true: done = false;
    return done;
  }
};

/*
 Combine termination conditions
 */

enum HowTermCombine { ALL, ONE };

class CombineTerminate : public BaseTerminate
{
protected:
  vector<BaseTerminate> terms_;
  HowTermCombine howCombine_;
  
public:
  CombineTerminate(vector<BaseTerminate> terms, HowTermCombine how_combine)
  :BaseTerminate(), howCombine_(how_combine)
  {
  }
  
  const bool is_terminated(shared_ptr<System> _sys) const
  {
    bool done;
    howCombine_== ALL ? done=true: done=false;
    for (int i = 0; i < terms_.size(); i++)
    {
      if (terms_[i].is_terminated(_sys) == ! done)
      {
        done=!done;
        break;
      }
    }
    return done;
  }
  
};

/*
 Class for implementing Brownian Dynamics on a a system
 */
class BD
{
protected:
  vector<double> transDiffConsts_;  // translational diffusion constants
  vector<double> rotDiffConsts_;  // rotational diffusion constants
  
  bool diff_; // include random kicks in dynamics
  bool force_; // include force calcs in dynamics
  double dt_;
  
  // random number generator object:
  mt19937 randGen_;
  shared_ptr<System> _sys_;
  shared_ptr<Constants> _consts_;
  
  // check if a molecule's new point causes it to collide with any other
  bool check_for_collision(int mol, Pt new_pt);
  
  // updates on individual molecules:
  void indi_trans_update(int i, MyVector<double> fi);
  void indi_rot_update(int i, MyVector<double> tau_i);
  
  // compute timestep for BD
  double compute_dt( double dist );
  
  // compute the smallest distance between two molecule centers
  double compute_min_dist( );
  
  // return a random vector with each element drawn from a Gaussian
  Pt rand_vec(double mean, double var);
  
  // update System time
  void update_sys_time(double dt) { _sys_->set_time(_sys_->get_time() + dt); }
  
public:
  BD(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
     vector<double> trans_diff_consts,
     vector<double> rot_diff_consts,
     bool diff = true, bool force = true);
  
  // update the system with Brownian dynamics given forces and torques on every
  // molecule
  void bd_update(VecOfVecs<double>::type F, VecOfVecs<double>::type tau);
  
  System get_system() { return *_sys_; }
  
};



#endif /* BD_h */

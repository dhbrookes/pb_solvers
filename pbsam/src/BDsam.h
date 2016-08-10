//
//  BDsam.h
//  pbsam_xcode
//
//  Created by David Brookes on 8/8/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BDsam_h
#define BDsam_h

#include <stdio.h>
#include <vector>
#include <random>
#include <memory>
#include "System.h"

using namespace std;
/*
 Class for performing a brownian dynamics step
 */
class BDStep
{
protected:
  vector<double> transDiffConsts_;  // translational diffusion constants
  vector<double> rotDiffConsts_;  // rotational diffusion constants
  
  bool diff_; // include random kicks in dynamics
  bool force_; // include force calcs in dynamics
  double dt_;
  double min_dist_;
  
  // random number generator object:
  mt19937 randGen_;
  shared_ptr<System> _sys_;
  shared_ptr<Constants> _consts_;
  
  // check if a molecule's new point causes it to collide with any other
  bool check_for_collision(int mol, Pt new_pt);
  
  // updates on individual molecules:
  void indi_trans_update(int i, Pt fi);
  void indi_rot_update(int i, Pt tau_i);
  
  // compute timestep for BD
  double compute_dt( );
  
  // compute the smallest distance between two molecule centers
  void compute_min_dist( );
  
  // return a random vector with each element drawn from a Gaussian
  Pt rand_vec(double mean, double var);
  
  // update System time
  void update_sys_time(double dt) { _sys_->set_time(_sys_->get_time() + dt); }
  
public:
  BDStep(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
         vector<double> trans_diff_consts,
         vector<double> rot_diff_consts,
         bool diff = true, bool force = true);
  
  // Constructor where diffusion constants are read from system:
  BDStep(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
         bool diff = true, bool force = true);
  
  // update the system with Brownian dynamics given forces and torques on every
  // molecule
  void bd_update(shared_ptr<vector<Pt> > _F,
                 shared_ptr<vector<Pt> > _tau);
  
  shared_ptr<System> get_system() { return _sys_; }
  double get_dt()                 { return dt_; }
  double get_min_dist()           { return min_dist_; }
  
};

#endif /* BDsam_h */

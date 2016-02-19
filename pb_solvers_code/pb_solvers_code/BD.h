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
  
  // random number generator object:
  mt19937 randGen_;
  shared_ptr<System> _sys_;
  
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
  
public:
  BD(System sys, vector<double> trans_diff_consts,
     vector<double> rot_diff_consts, bool diff = true, bool force = true);
  
  // update the system with Brownian dynamics given forces and torques on every
  // molecule
  void bd_update(VecOfVecs<double>::type F, VecOfVecs<double>::type tau);
  
  System get_system() { return *_sys_; }
  
};


#endif /* BD_h */

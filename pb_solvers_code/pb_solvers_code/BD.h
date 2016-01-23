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
  double dt_;  // time step
  vector<double> transDiffConsts_;  // translational diffusion constants
  vector<double> rotDiffConsts_;  // rotational diffusion constants
  
  // random number generator object:
  mt19937 randGen_;
  
  shared_ptr<System> _sys_;
  
  double kb_;
  double T_;  // temperature
  
  // check if a molecule's new point causes it to collide with any other
  bool check_for_collision(int mol, Pt new_pt);
  
  // updates on individual molecules:
  void indi_trans_update(int i, MyVector<double> fi);
  void indi_rot_update(int i, MyVector<double> tau_i);
  
  // return a random vector with each element drawn from a Gaussian
  Pt rand_vec(double mean, double var);
  
public:
  BD(System sys, double dt, vector<double> trans_diff_consts,
     vector<double> rot_diff_consts, double kb, double T);
  
  // update the system with Brownian dynamics given forces and torques on every
  // molecule
  void bd_update(VecOfVecs<double>::type F,
                 VecOfVecs<double>::type tau);
  
};


#endif /* BD_h */

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
#include <random>
#include <memory>
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
 Class for contact based termination. This terminates based on whether
 the specified molecule molecule pair is within a given cutoff of each other
 */
class ContactTerminate : public BaseTerminate
{
protected:
  double dist_contact_; //termination time
  int mol1_;
  int mol2_;
  
public:
  ContactTerminate(vector<int> mol, double distance)
  :BaseTerminate(), mol1_(mol[0]), mol2_(mol[1]), dist_contact_(distance)
  {
  }
  
  const bool is_terminated(shared_ptr<System> _sys) const
  {
    int i, j, idx1, idx2;
    for ( i = 0; i < _sys->get_typect(mol1_); i++)
    {
      for ( j = 0; j < _sys->get_typect(mol2_); j++)
      {
        idx1 = _sys->get_mol_global_idx( mol1_, i);
        idx2 = _sys->get_mol_global_idx( mol2_, j);

        double mol_c2c = _sys->get_pbc_dist_vec(idx1, idx2).norm();
    
        if (mol_c2c <= dist_contact_) return true;
      }
    }
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
enum BoundaryType { LEQ, GEQ };

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
    bool done = false;
    int i;
    for ( i = 0; i < _sys->get_typect(molIdx_); i++)
    {
      Pt mol_coord = _sys->get_centeri(_sys->get_mol_global_idx( molIdx_, i));
      double test_val;
      if (coordType_ == X)      test_val = mol_coord.x();
      else if (coordType_ == Y) test_val = mol_coord.y();
      else if (coordType_ == Z) test_val = mol_coord.z();
      else                      test_val = mol_coord.norm();
      
      if (boundType_ == LEQ)
        test_val <= boundaryVal_ ? done = true: done = false;
      else if (boundType_ == GEQ)
        test_val >= boundaryVal_ ? done = true: done = false;
    }
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


/*
 Class for running a full BD simulation
 */
class BDRun
{
protected:
  shared_ptr<BDStep> _stepper_;
  shared_ptr<ASolver> _asolver_;
  shared_ptr<BasePhysCalc> _physCalc_;
  shared_ptr<BaseTerminate> _terminator_;
  int maxIter_;
  double prec_;
  
public:
  // num is the number of bodies to perform calculations on (2, 3 or all).
  // If num=0, then the equations will be solved exactly
  BDRun(shared_ptr<ASolver> _asolv, shared_ptr<BaseTerminate> _terminator,
        int num=0, bool diff = true, bool force = true, int maxiter=1000,
        double prec=1e-4);
  
  void run();
};




#endif /* BD_h */

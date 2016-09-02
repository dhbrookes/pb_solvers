//
//  BaseBD.hpp
//  pbam_xcode
//
//  Created by David Brookes on 9/1/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BaseBD_h
#define BaseBD_h

#include <stdio.h>
#include <random>
#include <memory>
#include <vector>
#include "BaseSys.h"
#include "BasePhysCalc.h"

using namespace std;
/*
 Base class for implementing termination conditions in BD
 */
class BaseTerminate
{
public:
  BaseTerminate() { }
  
  virtual const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  virtual string get_how_term(shared_ptr<BaseSystem> _sys);
};

enum CoordType { X, Y, Z, R };
enum BoundaryType { LEQ, GEQ };

/*
 Class for time based termination
 */
class TimeTerminate : public BaseTerminate
{
protected:
  double endTime_; //termination time
  string how_term_;
  
public:
  TimeTerminate(double end_time);
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
};

/*
 Class for coordinate based termination. This terminates based on whether
 the specified MoleculeAM satisfies the BoundaryType condition on the CoordType
 with the given boundary value.
 */
class CoordTerminate : public BaseTerminate
{
protected:
  double boundaryVal_;
  int molIdx_;
  CoordType coordType_;
  BoundaryType boundType_;
  string how_term_;
  
public:
  CoordTerminate(int mol_idx, CoordType coord_type,
                 BoundaryType bound_type, double bound_val);
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
};

/*
 Combine termination conditions
 */
enum HowTermCombine { ALL, ONE };

class CombineTerminate: public BaseTerminate
{
protected:
  vector<shared_ptr<BaseTerminate> > terms_;
  HowTermCombine howCombine_;
  
public:
  CombineTerminate(vector<shared_ptr<BaseTerminate> > terms,
                   HowTermCombine how_combine);
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  
  string get_how_term(shared_ptr<BaseSystem> _sys);
};


/*
 Base class for performing a brownian dynamics step
 */
class BaseBDStep
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
  shared_ptr<BaseSystem> _sys_;
  shared_ptr<Constants> _consts_;
  
  // check if a molecule's new point causes it to collide with any other
  //  bool check_for_collision(int mol, Pt new_pt);
  
  // updates on individual moleculess:
  void indi_trans_update(int i, Pt fi);
  void indi_rot_update(int i, Pt tau_i);
  
  // compute timestep for BD
  double compute_dt( );
  
  // compute the smallest distance between two molecule centers
  // this is really the only function that will differ in sub classes
  virtual void compute_min_dist( ) { }
  
  // return a random vector with each element drawn from a Gaussian
  Pt rand_vec(double mean, double var);
  
  // update System time
  void update_sys_time(double dt) { _sys_->set_time(_sys_->get_time() + dt); }
  
public:
  BaseBDStep(shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
           vector<double> trans_diff_consts,
           vector<double> rot_diff_consts,
           bool diff = true, bool force = true);
  
  // Constructor where diffusion constants are read from system:
  BaseBDStep (shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
           bool diff = true, bool force = true);
  
  // update the system with Brownian dynamics given forces and torques on every
  // molecule
  void bd_update(shared_ptr<vector<Pt> > _F,
                 shared_ptr<vector<Pt> > _tau);
  
  shared_ptr<BaseSystem> get_system() { return _sys_; }
  double get_dt()                     { return dt_; }
  double get_min_dist()               { return min_dist_; }
};

/*
 Class for running a full BD simulation
 */
class BaseBDRun
{
protected:
  shared_ptr<BaseBDStep> _stepper_;
  shared_ptr<BasePhysCalc> _physCalc_;
  shared_ptr<BaseTerminate> _terminator_;
  
  string outfname_; //outputfile
  
  int maxIter_;
  double prec_;
  
public:
  // num is the number of bodies to perform calculations on (2, 3 or all).
  // If num=0, then the equations will be solved exactly
  BaseBDRun(shared_ptr<BaseTerminate> _terminator,
          string outfname, int num=0, bool diff = true, bool force = true,
          int maxiter=1e8, double prec=1e-4);
  
  virtual void run(string xyzfile = "test.xyz", string statfile = "stats.dat",
                   int nSCF = 0) { }
  
  Pt get_force_i(int i)      {return _physCalc_->get_forcei(i);}
  Pt get_torque_i(int i)     {return _physCalc_->get_taui(i);}
  double get_energy_i(int i) {return _physCalc_->calc_ei(i);}
  
  Pt get_forcei_conv(int i)      {return _physCalc_->get_forcei(i);}
  Pt get_torquei_conv(int i)     {return _physCalc_->get_taui(i);}
  double get_energyi_conv(int i) {return _physCalc_->calc_ei(i);}
};




#endif /* BaseBD_h */

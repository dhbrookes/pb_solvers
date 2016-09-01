//
//  BDAM.h
//  pb_solvers_code
//
/*
Copyright (c) 2015, Teresa Head-Gordon, Lisa Felberg, Enghui Yap, David Brookes
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
* Neither the name of UC Berkeley nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef BDAM_h
#define BDAM_h

#include <stdio.h>
#include <random>
#include <memory>
#include "ElectrostaticsAM.h"
#include "BaseBD.h"


/*
 Class for contact based termination. This terminates based on whether
 the specified molecule-molecule pair is within a given cutoff of each other
 */
class ContactTerminateAM : public BaseTerminate
{
protected:
  double dist_contact_; //termination time
  int mol1_;
  int mol2_;
  string how_term_;
  
public:
  ContactTerminateAM(vector<int> mol, double distance);
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
};

// a more complicated contact scheme based on distance between projections of
// atom positions on outer sphere
class ContactTerminateAM2 : public BaseTerminate
{
protected:
  int mol1_;
  int mol2_;
  
  double pad_; // distance between spheres if contact distance cannot be met
  
  vector<vector<int> > atPairs_;  // vector of size two vectors (atom index from each MoleculeAM type)
  vector<double> dists_;  // min distance between the above pairs
  
  string how_term_;
  
public:
  ContactTerminateAM2(vector<int> mol, vector<vector<int> > atpairs,
                      vector<double> dists, double pad);
  
  ContactTerminateAM2(ContactFile confile, double pad);
  
  void string_create();
  
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
};

/*
 Class for performing a brownian dynamics step
 */
class BDStepAM
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
  shared_ptr<SystemAM> _sys_;
  shared_ptr<Constants> _consts_;
  
  // check if a MoleculeAM's new point causes it to collide with any other
//  bool check_for_collision(int mol, Pt new_pt);
  
  // updates on individual MoleculeAMs:
  void indi_trans_update(int i, Pt fi);
  void indi_rot_update(int i, Pt tau_i);
  
  // compute timestep for BD
  double compute_dt( );
  
  // compute the smallest distance between two MoleculeAM centers
  void compute_min_dist( );
  
  // return a random vector with each element drawn from a Gaussian
  Pt rand_vec(double mean, double var);
  
  // update System time
  void update_sys_time(double dt) { _sys_->set_time(_sys_->get_time() + dt); }
  
public:
  BDStepAM(shared_ptr<SystemAM> _sys, shared_ptr<Constants> _consts,
     vector<double> trans_diff_consts,
     vector<double> rot_diff_consts,
     bool diff = true, bool force = true);
  
  // Constructor where diffusion constants are read from system:
  BDStepAM(shared_ptr<SystemAM> _sys, shared_ptr<Constants> _consts,
         bool diff = true, bool force = true);
  
  // update the system with Brownian dynamics given forces and torques on every
  // MoleculeAM
  void bd_update(shared_ptr<vector<Pt> > _F,
                 shared_ptr<vector<Pt> > _tau);
  
  shared_ptr<SystemAM> get_system() { return _sys_; }
  double get_dt()                 { return dt_; }
  double get_min_dist()           { return min_dist_; }
  
};


/*
 Class for running a full BD simulation
 */
class BDRunAM
{
protected:
  shared_ptr<BDStepAM> _stepper_;
  shared_ptr<ASolver> _asolver_;
  shared_ptr<BasePhysCalc> _physCalc_;
  shared_ptr<BaseTerminate> _terminator_;
  
  string outfname_; //outputfile
  
  int maxIter_;
  double prec_;
  
public:
  // num is the number of bodies to perform calculations on (2, 3 or all).
  // If num=0, then the equations will be solved exactly
  BDRunAM(shared_ptr<ASolver> _asolv, shared_ptr<BaseTerminate> _terminator,
        string outfname, int num=0, bool diff = true, bool force = true,
        int maxiter=1e8, double prec=1e-4);
  
  void run(string xyzfile = "test.xyz", string statfile = "stats.dat", 
           int nSCF = 0);

  Pt get_force_i(int i)      {return _physCalc_->get_forcei(i);}
  Pt get_torque_i(int i)     {return _physCalc_->get_taui(i);}
  double get_energy_i(int i) {return _physCalc_->calc_ei(i);}

  Pt get_forcei_conv(int i)      {return _physCalc_->get_forcei(i);}
  Pt get_torquei_conv(int i)     {return _physCalc_->get_taui(i);}
  double get_energyi_conv(int i) {return _physCalc_->calc_ei(i);}
};




#endif /* BDAM_h */

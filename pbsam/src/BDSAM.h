//
//  BDSAM.h
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

#ifndef BDSAM_h
#define BDSAM_h

#include <stdio.h>
#include <random>
#include <memory>
#include "BaseBD.h"
#include "ElectrostaticsSAM.h"


class ContactTerminateSAM : public BaseTerminate
{
protected:
  int mol1_;
  int mol2_;
  
  double pad_; // distance between spheres if contact distance cannot be met
  
  vector<vector<int> > atPairs_;  // vector of size two vectors (atom index from each MoleculeAM type)
  vector<double> dists_;  // min distance between the above pairs
  
  string how_term_;
  
public:
  ContactTerminateSAM(vector<int> mol, vector<vector<int> > atpairs,
                      vector<double> dists, double pad);
  
  ContactTerminateSAM(ContactFile confile, double pad);
  
  void string_create();
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
};

/*
 Class for performing a brownian dynamics step
 */
class BDStepSAM: public BaseBDStep
{
protected:

  shared_ptr<BaseSystem> _sys_;
  
  // compute the smallest distance between two MoleculeAM centers
  void compute_min_dist( );
  
public:
  BDStepSAM(shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
     vector<double> trans_diff_consts,
     vector<double> rot_diff_consts,
     bool diff = true, bool force = true);
  
  // Constructor where diffusion constants are read from system:
  BDStepSAM(shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
         bool diff = true, bool force = true);
  
};


/*
 Class for running a full BD simulation
 */
class BDRunSAM : public BaseBDRun
{
protected:
  shared_ptr<Solver> _solver_;
  shared_ptr<GradSolver> _gradSolv_;
  
public:
  // num is the number of bodies to perform calculations on (2, 3 or all).
  // If num=0, then the equations will be solved exactly
  BDRunSAM(shared_ptr<Solver> _solv, shared_ptr<GradSolver> _gradSolv,
        shared_ptr<BaseTerminate> _terminator, string outfname,
        bool diff = true, bool force = true,
        int maxiter = 2, double prec = 1e-4);
  
  void run(string xyzfile = "test.xyz", string statfile = "stats.dat", 
           int nSCF = 0);
};




#endif /* BD_h */

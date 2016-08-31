//
// PBAM.h
// pbam
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

#ifndef PBAM_H
#define PBAM_H

#include <memory>
#include <time.h>
#include "PBAMStruct.h"
#include "BD.h"


using namespace std;

class PBAM : protected PBAMInput
{
protected:
  shared_ptr<Setup> setp_;
  shared_ptr<System> syst_;
  shared_ptr<Constants> consts_;

  shared_ptr<BesselConstants> _bessl_consts_;
  shared_ptr<BesselCalc> _bessl_calc_;
  shared_ptr<SHCalcConstants> _sh_consts_;
  shared_ptr<SHCalc> _sh_calc_;

  double force_[MOL_MAX][3];
  double torque_[MOL_MAX][3];
  double nrg_intera_[MOL_MAX];

  int poles_;
  double solveTol_;

public:

  // Constructors
  PBAM();
  PBAM(string infile);
  // For APBS
  PBAM(const PBAMInput& pbami, vector<MoleculeAM> mls );

  friend PBAMInput getPBAMParams();

  // Copy constructors
  PBAM( const PBAM& pbam ) ;
  PBAM( const PBAM* pbam ) ;

  void check_setup();
  void check_system();

  void init_write_system();
  void initialize_coeff_consts();

  int run();
  // for running the APBS version
  PBAMOutput run_apbs( );

  void run_bodyapprox();
  void run_dynamics();
  void run_electrostatics();
  void run_energyforce();
};


#endif

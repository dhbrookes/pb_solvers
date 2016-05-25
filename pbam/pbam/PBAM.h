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
#include "setup.h"
#include "BD.h"
#include "Electrostatics.h"
#include <time.h>

using namespace std;

class PBAM
{
protected:
  shared_ptr<Setup> setp;
  shared_ptr<System> syst;
  shared_ptr<Constants> consts;

  int poles;
  double solveTol;
  
public:
  
  // Constructors
  PBAM();
  PBAM(string infile);
  PBAM(const struct PBAMInput& pbami ); // For APBS
  
  // for running the APBS version
  // struct PBAMOutput run( const AtomList& atomList);

  void get_check_inputs();

  int run();

  void run_bodyapprox();
  void run_dynamics();
  void run_electrostatics();
  void run_energyforce();
};


#endif

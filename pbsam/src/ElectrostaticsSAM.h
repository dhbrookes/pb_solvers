//
//  Electrostatics.h
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
#ifndef Electrostatics_h
#define Electrostatics_h

#include "PhysCalcSAM.h"
#include "BaseElectro.h"
#include <time.h>

#ifdef __OMP
#include <omp.h>
#endif

namespace pbsolvers
{

/*
 Class for printing out electrostatics of system
 */
class ElectrostaticSAM : public BaseElectro
{
protected:

  double eps_s;
  vector<shared_ptr<HMatrix> > _H_;

  double compute_pot_at( Pt point );

public:
  ElectrostaticSAM(vector<shared_ptr<HMatrix> > H, shared_ptr<SystemSAM> _sys,
                shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
                shared_ptr<Constants> _consts,
                int p, int npts = 150);

  ElectrostaticSAM(shared_ptr<Solver> solve, int npts=150);
};

} /* namespace pbsolvers */
#endif /* Electrostatics_h */

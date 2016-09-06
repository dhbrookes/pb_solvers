//
//  Electrostatics.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "ElectrostaticsSAM.h"

ElectrostaticSAM::ElectrostaticSAM(vector<shared_ptr<HMatrix> > H,
                             shared_ptr<SystemSAM> _sys,
                             shared_ptr<SHCalc> _shCalc,
                             shared_ptr<BesselCalc> _bCalc,
                             shared_ptr<Constants> _consts,
                             int p, int npts)
: BaseElectro(_sys, _shCalc, _bCalc, _consts, p, npts), _H_(H)
{
  compute_pot();
}

ElectrostaticSAM::ElectrostaticSAM(shared_ptr<Solver> solve, int npts)

:BaseElectro(solve->get_sys(), solve->get_sh(), solve->get_bessel(),
             solve->get_consts(), solve->get_p(), npts),
_H_(solve->get_all_H())
{
  compute_pot();
}

double ElectrostaticSAM::compute_pot_at( Pt point )
{

  int mol, sph, Nmol = _sys_->get_n();
  double rad, pot    = 0.0;
  Pt center, dist;
  MyMatrix<cmplx> localK;
  
  for ( mol = 0; mol < Nmol; mol++)
  {
    for ( sph = 0; sph < _sys_->get_Ns_i(mol); sph++)
    {
      center = _sys_->get_centerik(mol, sph);
      rad    = _sys_->get_aik(mol, sph);
      dist   = point - center;
      localK = get_local_exp(dist, rad);
      pot += lotan_inner_prod( _H_[mol]->get_mat_k(sph), localK, p_);
    }
  }
  return pot;
}


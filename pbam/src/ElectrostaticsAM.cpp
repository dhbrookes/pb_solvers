//
//  ElectrostaticsAM.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "ElectrostaticsAM.h"

namespace pbsolvers
{

ElectrostaticAM::ElectrostaticAM(shared_ptr<VecOfMats<cmplx>::type> _A,
                             shared_ptr<SystemAM> _sys,
                             shared_ptr<SHCalc> _shCalc,
                             shared_ptr<BesselCalc> _bCalc,
                             shared_ptr<Constants> _consts,
                             int p, int npts)
: BaseElectro(_sys, _shCalc, _bCalc, _consts, p, npts), _A_(_A)
{
  compute_pot();
}

ElectrostaticAM::ElectrostaticAM(shared_ptr<ASolver> solve, int npts)
: BaseElectro(solve->get_sys(), solve->get_sh(), solve->get_bessel(),
              solve->get_consts(), solve->get_p(), npts),
_A_(solve->get_A())
{
  compute_pot();
}

double ElectrostaticAM::compute_pot_at( Pt point )
{

  int mol, Nmol      = _sys_->get_n();
  double rad, pot    = 0.0;
  Pt center, dist;
  MyMatrix<cmplx> localK;
  
  for ( mol = 0; mol < Nmol; mol++)
  {
    center = _sys_->get_centerik(mol, 0);
    rad    = _sys_->get_aik(mol, 0);
    dist   = point - center;
    localK = get_local_exp(dist, _sys_->get_lambda());
    pot += lotan_inner_prod( _A_->operator[](mol), localK, p_);
  }
 
  return pot;
}

} /* namespace pbsolvers */

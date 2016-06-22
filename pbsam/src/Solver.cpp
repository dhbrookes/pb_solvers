//
//  Solver.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solver.h"


Solver::Solver(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
               shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
               int p)
:_E_(_sys->get_n()),
_LE_(_sys->get_n()),
_IE_(_sys->get_n()),
_LF_(_sys->get_n()),
_LH_(_sys->get_n()),
_LHN_(_sys->get_n()),
_XF_(_sys->get_n()),
_XH_(_sys->get_n()),
_H_(_sys->get_n()),
_F_(_sys->get_n()),
_prevH_(_sys->get_n()),
_prevF_(_sys->get_n()),
_shCalc_(_shCalc),
_bCalc_(_bCalc),
_consts_(_consts),
_sys_(_sys),
p_(p),
kappa_(_consts->get_kappa())
{
  _reExConsts_ = make_shared<ReExpCoeffsConstants> (kappa_,
                                                    _sys_->get_lambda(), p_);
  _T_ = make_shared<TMatrix> (p_, _sys_, _shCalc_, _consts_,
                             _bCalc_, _reExConsts_);
  
  _expConsts_ = make_shared<ExpansionConstants>(p_);
  shared_ptr<Molecule> _mol;
  // intialize all matrices
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    _mol = make_shared<Molecule>(_sys_->get_molecule(I));
    double kappa = 0.0;
    
    _E_[I] = make_shared<EMatrix> (I, _sys_->get_Ns_i(I), p_);
    _E_[I]->calc_vals((*_mol), _shCalc_, _consts_->get_dielectric_prot());
    
    _LE_[I] = make_shared<LEMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LE_[I]->calc_vals((*_mol), _shCalc_, _consts_->get_dielectric_prot());
    
    _H_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _H_[I]->init((*_mol), _shCalc_, _consts_->get_dielectric_prot());
    
    _F_[I] = make_shared<FMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    
    _prevH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _prevF_[I] = make_shared<FMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    
    _IE_[I] = make_shared<IEMatrix>(I, _mol, _shCalc, p_, _expConsts_, false);
    _IE_[I]->calc_vals(_mol, _shCalc_);
    
    _LF_[I] = make_shared<LFMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LH_[I] = make_shared<LHMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
    _LHN_[I] = make_shared<LHNMatrix> (I, _sys_->get_Ns_i(I), p_);
    
    _XF_[I] = make_shared<XFMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _consts_->get_dielectric_prot(),
                                     _consts_->get_dielectric_water(),
                                     (*_mol), _E_[I], _LE_[I]);
    _XH_[I] = make_shared<XHMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     (*_mol), _E_[I], _LE_[I]);
  
    update_prev();
  }
}

double Solver::iter()
{
  double mu = 0;
  Molecule mol;
  // start an iteration
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    mol = _sys_->get_molecule(I);
    
    _LH_[I]->init(_sys_->get_molecule(I),_H_[I],_shCalc_,_bCalc_,_expConsts_);
    _LH_[I]->calc_vals(_T_, _H_[I], _shCalc_, _sys_, _bCalc_);
    
    _LF_[I]->calc_vals(_T_, _F_[I], _shCalc_, _sys_);
    
    _XF_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_,
                      _LH_[I], _LF_[I], _LHN_[I], 0.0);
    _XH_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_, _LH_[I],
                       _LF_[I], _LHN_[I], 0.0);
    
    _F_[I]->calc_vals(mol, _prevF_[I], _XF_[I], _prevH_[I], _IE_[I], _bCalc_);
    _H_[I]->calc_vals(mol, _prevH_[I], _XH_[I], _prevF_[I], _IE_[I], _bCalc_);
    
    mu += HMatrix::calc_converge(_H_[I], _prevH_[I]);
  }
  
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    _LHN_[I]->calc_vals(_T_, _H_); // requires all Hs so do it after they are
                                    // calculated
  }
  update_prev();
  return mu;
}


void Solver::update_prev()
{
  for (int I = 0; I < _H_.size(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      for (int n = 0; n < p_; n++)
      {
        for (int m = 0; m < p_; m++)
        {
          _prevH_[I]->set_mat_knm(k, n, m, _H_[I]->get_mat_knm(k, n, m));
          _prevF_[I]->set_mat_knm(k, n, m, _F_[I]->get_mat_knm(k, n, m));
        }
      }
    }
  }
}

void Solver::solve(double tol, int maxiter)
{
  double mu;
  for (int t = 0; t < maxiter; t++)
  {
    mu = iter();
    if (mu < tol) break;
  }
}


void Solver::solve_inner()
{
  double mu;
  double solveTol = 1e-6;
  int maxiter = 20;
  for (int t = 0; t < maxiter; t++)
  {
    mu = iter();
    if (mu < solveTol) break;
  }
}


void Solver::reset_all()
{
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    _E_[I]->reset_mat();
    _LE_[I]->reset_mat();
    _IE_[I]->reset_mat();
    _LF_[I]->reset_mat();
    _LH_[I]->reset_mat();
    _LHN_[I]->reset_mat();
    _XF_[I]->reset_mat();
    _XH_[I]->reset_mat();
    _H_[I]->reset_mat();
    _F_[I]->reset_mat();
    _prevH_[I]->reset_mat();
    _prevF_[I]->reset_mat();
  }
}


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
               shared_ptr<TMatrix> _T, int p)
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
_shCalc_(_shCalc),
_bCalc_(_bCalc),
_consts_(_consts),
_sys_(_sys),
p_(p),
kappa_(_consts->get_kappa()),
_T_(_T)
{
//  _reExConsts_ = make_shared<ReExpCoeffsConstants> (kappa_,
//                                                    _sys_->get_lambda(), p_);
//  _T_ = make_shared<TMatrix> (p_, _sys_, _shCalc_, _consts_,
//                             _bCalc_, _reExConsts_);
  
  auto _expConst = make_shared<ExpansionConstants>(p_);
  shared_ptr<Molecule> _mol;
  // intialize all matrices
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    _mol = make_shared<Molecule> (_sys_->get_molecule(I));
    
    _E_[I] = make_shared<EMatrix> (I, _sys_->get_Ns_i(I), p_);
    _E_[I]->calc_vals((*_mol), _shCalc_, _consts_->get_dielectric_prot());
    
    _LE_[I] = make_shared<LEMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LE_[I]->calc_vals((*_mol), _shCalc_, _consts_->get_dielectric_prot());
    
    _IE_[I] = make_shared<IEMatrix> (I, _sys_->get_Ns_i(I), p_, _expConst);
    _IE_[I]->calc_vals(_mol, _shCalc_);
    
    _LF_[I] = make_shared<LFMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LH_[I] = make_shared<LHMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
    _LHN_[I] = make_shared<LHNMatrix> (I, _sys_->get_Ns_i(I), p_);
    
    _XF_[I] = make_shared<XFMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _consts_->get_dielectric_prot(),
                                     _consts_->get_dielectric_water());
    _XH_[I] = make_shared<XHMatrix> (I, _sys_->get_Ns_i(I), p_);
    
    _H_[I] = make_shared<HMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
    _F_[I] = make_shared<FMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
    _prevH_[I] = make_shared<HMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
    _prevF_[I] = make_shared<FMatrix> (I, _sys_->get_Ns_i(I), p_, kappa_);
  }
}

double Solver::iter()
{
  double mu = 0;
  Molecule mol;
  // start an iteration by calculating XF and XH (is this right?)
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    mol = _sys_->get_molecule(I);
    _XF_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_, _E_[I], _LE_[I],
                      _LH_[I], _LF_[I], _LHN_[I], kappa_);
    _XH_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_, _E_[I], _LE_[I],
                      _LH_[I], _LF_[I], _LHN_[I]);
    
    _H_[I]->calc_vals(mol, _prevH_[I], _XH_[I], _prevF_[I], _IE_[I], _bCalc_);
    _F_[I]->calc_vals(mol, _prevF_[I], _XF_[I], _prevH_[I], _IE_[I], _bCalc_);
    
    _LF_[I]->calc_vals(_T_, _F_[I], _shCalc_, _sys_);
    _LH_[I]->calc_vals(_T_, _H_[I], _shCalc_, _sys_, _bCalc_);
    
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



GradSolver::GradSolver(shared_ptr<System> _sys,
                       shared_ptr<Constants> _consts,
                       shared_ptr<SHCalc> _shCalc,
                       shared_ptr<PrecalcBessel> _bCalc,
                       shared_ptr<TMatrix> _T,
                       vector<shared_ptr<FMatrix> > _F,
                       vector<shared_ptr<HMatrix> > _H,
                       int p)
:_F_(_F), _H_(_H), _T_(_T), _bCalc_(_bCalc), _shCalc_(_shCalc),
_sys_(_sys), _consts_(_consts), kappa_(_consts->get_kappa()),
dF_(_sys->get_n(), vector<shared_ptr<GradFMatrix> > (_sys->get_n())),
dH_(_sys->get_n(), vector<shared_ptr<GradHMatrix> > (_sys->get_n())),
dWF_(_sys->get_n(), vector<shared_ptr<GradWFMatrix> > (_sys->get_n())),
dWH_(_sys->get_n(), vector<shared_ptr<GradWHMatrix> > (_sys->get_n())),
dLF_(_sys->get_n(), vector<shared_ptr<GradLFMatrix> > (_sys->get_n())),
dLH_(_sys->get_n(), vector<shared_ptr<GradLHMatrix> > (_sys->get_n())),
dLHN_(_sys->get_n(), vector<shared_ptr<GradLHNMatrix> > (_sys->get_n()))
{
  dF_.reserve(_sys_->get_n());
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    for (int J = 0; J < _sys_->get_n(); J++)
    {
      dF_[I][J] = make_shared<GradFMatrix> (I, J, _sys_->get_Ns_i(I), p_);
      dWF_[I][J] = make_shared<GradWFMatrix> (I, J, _sys_->get_Ns_i(I), p_,
                                              _consts->get_dielectric_prot(),
                                              _consts->get_dielectric_water(),
                                              _consts->get_kappa());
      
      dLF_[I][J] = make_shared<GradLFMatrix> (I, J, _sys_->get_Ns_i(I), p_);
      dLHN_[I][J] = make_shared<GradLHNMatrix> (I, J, _sys_->get_Ns_i(I), p_);
      
      dH_[I][J] = make_shared<GradHMatrix> (I,J,_sys_->get_Ns_i(I),p_,kappa_);
      dWH_[I][J] = make_shared<GradWHMatrix> (I,J,_sys_->get_Ns_i(I),p_,kappa_);
      dLH_[I][J] = make_shared<GradLHMatrix> (I,J,_sys_->get_Ns_i(I),p_,kappa_);
      
    }
  }
  
}



void GradSolver::solve()
{
  shared_ptr<Molecule> molI;
  for (int J = 0; J < _sys_->get_n(); J++)  // with respect to
  {
    for (int I = 0; I < _sys_->get_n(); I++)
    {
      molI = make_shared<Molecule>(_sys_->get_molecule(I));
      for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      {
        dLF_[J][I]->calc_val_k(k, molI, _shCalc_, _T_, dF_[J][I]);
        dLH_[J][I]->calc_val_k(k, molI, _bCalc_, _shCalc_, _T_, dH_[J][I]);
        dLHN_[J][I]->calc_val_k(k, _sys_, _T_, _H_, dH_[J]);

        dWF_[J][I]->calc_val_k(k, molI, _bCalc_, dH_[J][I], dF_[J][I],
                               dLH_[J][I], dLHN_[J][I], dLF_[J][I]);
        dWH_[J][I]->calc_val_k(k, molI, _bCalc_, dH_[J][I], dF_[J][I],
                               dLH_[J][I], dLHN_[J][I], dLF_[J][I]);
        
        dF_[J][I]->calc_val_k(k, _IE_[I], dWF_[J][I]);
        dH_[J][I]->calc_val_k(k, _bCalc_, _IE_[I], dWH_[J][I]);
      }
    }
  }
}


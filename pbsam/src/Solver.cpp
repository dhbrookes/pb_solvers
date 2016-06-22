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
    _mol = _sys_->get_molecule(I);
    double kappa = 0.0;
    
    _E_[I] = make_shared<EMatrix> (I, _sys_->get_Ns_i(I), p_);
    _E_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());
    
    _LE_[I] = make_shared<LEMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LE_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());
    
    _H_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _H_[I]->init(_mol, _shCalc_, _consts_->get_dielectric_prot());
    
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
                                     _mol, _E_[I], _LE_[I]);
    _XH_[I] = make_shared<XHMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _mol, _E_[I], _LE_[I]);
    
    update_prev_all();
  }
}

double Solver::calc_converge_H(int I, int k, bool inner)
{
  double mu(0), rt(0), num(0), den(0);
  cmplx hnm_curr, hnm_prev;
  
  for (int n = 0; n < _H_[I]->get_p(); n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      hnm_curr = _H_[I]->get_mat_knm(k, n, m);
      hnm_prev = _prevH_[I]->get_mat_knm(k, n, m);
      
      if (inner)
      {
        if ( norm(hnm_prev) < 1e-15 ) // checking for mat values = 0
        {
          if ( norm(hnm_curr) < 1e-15 )
            continue;
          else
            rt = 1.0;
        } else
        {
          if ( norm(hnm_curr) < 1e-15 )
            rt = 1.0;
          else
            rt = norm(hnm_curr - hnm_prev)/(norm(hnm_curr) + norm(hnm_prev));
        }
        mu += rt;
      }
      else
      {
        num += norm(hnm_curr - hnm_prev);
        den += norm(hnm_curr) + norm(hnm_prev);
      }
    }
  }
  
  if (inner)
    return mu / double(4.0*p_*p_);
  else
  {
    if ( fabs(den) > 1e-15 )
      mu = (num / den);
    else
      mu = 0.0;
    return mu / 4.0;
  }
    
    
}

double Solver::iter()
{
  double devk, mu = 0;
  shared_ptr<Molecule> mol;
  // start an iteration
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    mol = _sys_->get_molecule(I);
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      double dev = 1e20;
      double tol = 1e-6;
      int ct(0), mxCt(20);
      _LH_[I]->init(_sys_->get_molecule(I),_H_[I],_shCalc_,_bCalc_,_expConsts_);
      _LH_[I]->calc_vals(_T_, _H_[I], k);
      
      _LF_[I]->calc_vals(_T_, _F_[I], _shCalc_, _sys_, k);
      
      _XF_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_,
                         _LH_[I], _LF_[I], _LHN_[I], 0.0, k);
      _XH_[I]->calc_vals(_sys_->get_molecule(I), _bCalc_, _LH_[I],
                         _LF_[I], _LHN_[I], 0.0, k);
      
//      _F_[I]->calc_vals(mol,_prevF_[I],_XF_[I],_prevH_[I],_IE_[I],_bCalc_,k);
//      _H_[I]->calc_vals(mol,_prevH_[I],_XH_[I],_F_[I],_IE_[I],_bCalc_,k);
//      
//      dev = calc_converge_H(I,k);
      
      while ( dev > tol && ct < mxCt )
      {
        _F_[I]->calc_vals(mol,_F_[I],_XF_[I],_prevH_[I],_IE_[I],_bCalc_,k);
        _H_[I]->calc_vals(mol,_prevH_[I],_XH_[I],_F_[I],_IE_[I],_bCalc_,k);
        
        dev = calc_converge_H(I,k,true);
        update_prevH(I, k);
        
        ct++;
      }
      
      cout << "This is molecule " << k << endl;
      cout << "This is F " ;
      _F_[I]->print_kmat(k);
      cout << "This is H ";
      _H_[I]->print_kmat(k);
      
      
      devk = calc_converge_H(I,k,false);
      if ( devk > mu )
      {
        cout << "This is mu " << mu << " and dev " << dev << endl;
        mu = dev;
//        cout << "This is curr " << endl;
//        _H_[I]->print_kmat(k);
//        cout << "This is prev" << endl;
//        _prevH_[I]->print_kmat(k);
        
      }
     
      for (int I = 0; I < _sys_->get_n(); I++)
      {
        for (int k = 0; k < _sys_->get_Ns_i(I); k++)
        {
          _LHN_[I]->calc_vals(_T_, _H_, k); // requires all Hs so do it at end
        }
      }
      
    }
  }

  return mu;
}


void Solver::update_prev_all()
{
  for (int I = 0; I < _H_.size(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      update_prevF(I, k);
      update_prevH(I, k);
    }
  }
}

void Solver::update_prevH(int I, int k)
{
  
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
      _prevH_[I]->set_mat_knm(k, n, m, _H_[I]->get_mat_knm(k, n, m));
//      _prevF_[I]->set_mat_knm(k, n, m, _F_[I]->get_mat_knm(k, n, m));
    }
  }
}

void Solver::update_prevF(int I, int k)
{
  
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
//      _prevH_[I]->set_mat_knm(k, n, m, _H_[I]->get_mat_knm(k, n, m));
        _prevF_[I]->set_mat_knm(k, n, m, _F_[I]->get_mat_knm(k, n, m));
    }
  }
}

void Solver::solve(double tol, int maxiter)
{
  double mu;
  for (int t = 0; t < maxiter; t++)
  {
    cout << "this is t " << t << endl;
    mu = iter();
    if (mu < tol) break;
  }
  
  cout << (*_H_[0]) << endl;
  
  cout << (*_F_[0]) << endl;
}


void Solver::reset_all()
{
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      _E_[I]->reset_mat(k);
      _LE_[I]->reset_mat(k);
      _LF_[I]->reset_mat(k);
      _LH_[I]->reset_mat(k);
      _LHN_[I]->reset_mat(k);
      _XF_[I]->reset_mat(k);
      _XH_[I]->reset_mat(k);
      _H_[I]->reset_mat(k);
      _F_[I]->reset_mat(k);
      _prevH_[I]->reset_mat(k);
      _prevF_[I]->reset_mat(k);
    }
    _IE_[I]->reset_mat();
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
      molI = _sys_->get_molecule(I);
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


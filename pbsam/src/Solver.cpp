//
//  Solver.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solver.h"

namespace pbsolvers
{

Solver::Solver(shared_ptr<SystemSAM> _sys, shared_ptr<Constants> _consts,
               shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
               int p, bool readImat, bool readHF,
               vector<vector<string> > imats,
               vector<vector<vector<string> > > expHF)
: Ns_tot_(0),
_E_(_sys->get_n()),
_LE_(_sys->get_n()),
_IE_(_sys->get_n()),
_LF_(_sys->get_n()),
_LH_(_sys->get_n()),
_LHN_(_sys->get_n()),
_XF_(_sys->get_n()),
_XH_(_sys->get_n()),
_H_(_sys->get_n()),
_F_(_sys->get_n()),
_outerH_(_sys->get_n()),
_rotH_(_sys->get_n()),
_prevH_(_sys->get_n()),
_shCalc_(_shCalc),
_bCalc_(_bCalc),
_consts_(_consts),
_sys_(_sys),
dev_sph_Ik_(_sys->get_n()),
p_(p),
kappa_(_consts->get_kappa())
{
  _reExConsts_ = make_shared<ReExpCoeffsConstants> (kappa_,
                                                    _sys_->get_lambda(), p_,
                                                    true);

  _precalcSH_ = make_shared<PreCalcSH>();

  _T_ = make_shared<TMatrix> (p_, _sys_, _shCalc_, _consts_,
                              _bCalc_, _reExConsts_);
  precalc_sh_lf_lh();
  precalc_sh_numeric();

  _expConsts_ = make_shared<ExpansionConstants>(p_);
  shared_ptr<BaseMolecule> _mol;
  // intialize all matrices
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    Ns_tot_ += _sys_->get_Ns_i(I);
    _mol = _sys_->get_moli(I);

    dev_sph_Ik_[I].resize(_mol->get_ns());
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      dev_sph_Ik_[I][k] = 1.0;

    double kappa = _consts->get_kappa();

    _E_[I] = make_shared<EMatrix> (I, _sys_->get_Ns_i(I), p_);
    _E_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());

    _LE_[I] = make_shared<LEMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LE_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());

    _H_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _F_[I] = make_shared<FMatrix>(I, _sys_->get_Ns_i(I), p_, 0.0);

    if (readHF)
    {
      for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      {
        _H_[I]->init_from_exp(expHF[I][k][0], k);
        _F_[I]->init_from_exp(expHF[I][k][1], k);
      }
    } else
    {
      _H_[I]->init(_mol, _shCalc_, _consts_->get_dielectric_prot());
    }
    _outerH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _prevH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _rotH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);

    _IE_[I] = make_shared<IEMatrix>(I, _mol, _shCalc, p_, _expConsts_, false);
    if (readImat)
    {
      for (int k = 0; k < _sys_->get_Ns_i(I); k++)
        _IE_[I]->init_from_file(imats[I][k], k);

    } else
      _IE_[I]->calc_vals(_mol, _shCalc_);

    _LF_[I] = make_shared<LFMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LH_[I] = make_shared<LHMatrix> (I, _sys_->get_Ns_i(I), p_, kappa);
    _LHN_[I] = make_shared<LHNMatrix> (I, _sys_->get_Ns_i(I), p_, _sys_);

    _XF_[I] = make_shared<XFMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _consts_->get_dielectric_prot(),
                                     _consts_->get_dielectric_water(),
                                     _mol, _E_[I], _LE_[I]);
    _XH_[I] = make_shared<XHMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _mol, _E_[I], _LE_[I]);
  }
  update_prev_all();
}


Solver::Solver(shared_ptr<SystemSAM> _sys, shared_ptr<Constants> _consts,
               shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
               int p, vector<shared_ptr<IEMatrix> > imats,
               vector<shared_ptr<HMatrix > > h_spol,
               vector<shared_ptr<FMatrix > > f_spol)
: Ns_tot_(0),
_E_(_sys->get_n()),
_LE_(_sys->get_n()),
_IE_(_sys->get_n()),
_LF_(_sys->get_n()),
_LH_(_sys->get_n()),
_LHN_(_sys->get_n()),
_XF_(_sys->get_n()),
_XH_(_sys->get_n()),
_H_(_sys->get_n()),
_F_(_sys->get_n()),
_outerH_(_sys->get_n()),
_rotH_(_sys->get_n()),
_prevH_(_sys->get_n()),
_shCalc_(_shCalc),
_bCalc_(_bCalc),
_consts_(_consts),
_sys_(_sys),
dev_sph_Ik_(_sys->get_n()),
p_(p),
kappa_(_consts->get_kappa())
{
  int molt;
  _reExConsts_ = make_shared<ReExpCoeffsConstants> (kappa_,
                                                    _sys_->get_lambda(),
                                                    p_, true);
  _precalcSH_ = make_shared<PreCalcSH>();

  _T_ = make_shared<TMatrix> (p_, _sys_, _shCalc_, _consts_,
                              _bCalc_, _reExConsts_);
  precalc_sh_lf_lh();
  precalc_sh_numeric();

  _expConsts_ = make_shared<ExpansionConstants>(p_);
  shared_ptr<BaseMolecule> _mol;
  // intialize all matrices
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    Ns_tot_ += _sys_->get_Ns_i(I);
    _mol = _sys_->get_moli(I);
    molt = _mol->get_type();

    dev_sph_Ik_[I].resize(_mol->get_ns());
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      dev_sph_Ik_[I][k] = 1.0;

    double kappa = _consts->get_kappa();

    _E_[I] = make_shared<EMatrix> (I, _sys_->get_Ns_i(I), p_);
    _E_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());

    _LE_[I] = make_shared<LEMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LE_[I]->calc_vals(_mol, _shCalc_, _consts_->get_dielectric_prot());

    _IE_[I] = make_shared<IEMatrix>(I, _mol, _shCalc, p_, _expConsts_, false);
    _IE_[I]->init_from_other(imats[I]);

    _H_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _F_[I] = make_shared<FMatrix>(I, _sys_->get_Ns_i(I), p_, 0.0);
    if (h_spol.size() == 0)
    {
      _H_[I]->init(_mol, _shCalc_, _consts_->get_dielectric_prot());
    } else
    {
      _H_[I]->set_all_mats(h_spol[molt]);
      _F_[I]->set_all_mats(f_spol[molt]);
    }

    _outerH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _prevH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);
    _rotH_[I] = make_shared<HMatrix>(I, _sys_->get_Ns_i(I), p_, kappa);

    _LF_[I] = make_shared<LFMatrix> (I, _sys_->get_Ns_i(I), p_);
    _LH_[I] = make_shared<LHMatrix> (I, _sys_->get_Ns_i(I), p_, kappa);
    _LHN_[I] = make_shared<LHNMatrix> (I, _sys_->get_Ns_i(I), p_, _sys_);

    _XF_[I] = make_shared<XFMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _consts_->get_dielectric_prot(),
                                     _consts_->get_dielectric_water(),
                                     _mol, _E_[I], _LE_[I]);
    _XH_[I] = make_shared<XHMatrix> (I, _sys_->get_Ns_i(I), p_,
                                     _mol, _E_[I], _LE_[I]);
  }
  update_prev_all();
}


Solver::Solver(shared_ptr<Solver> solvin)
:
p_(solvin->p_), kappa_(solvin->kappa_),
Ns_tot_(solvin->Ns_tot_),
_E_(solvin->_E_),
_LE_(solvin->_LE_),
_IE_(solvin->_IE_),
_LF_(solvin->_LF_),
_LH_(solvin->_LH_),
_LHN_(solvin->_LHN_),
_XF_(solvin->_XF_),
_XH_(solvin->_XH_),
_H_(solvin->_H_),
_outerH_(solvin->_outerH_),
_rotH_(solvin->_rotH_),
_prevH_(solvin->_prevH_),
_F_(solvin->_F_),
_T_(solvin->_T_),
_shCalc_(solvin->_shCalc_),
_bCalc_(solvin->_bCalc_),
_consts_(solvin->_consts_),
_sys_(solvin->_sys_),
_reExConsts_(solvin->_reExConsts_),
_expConsts_(solvin->_expConsts_),
dev_sph_Ik_(solvin->dev_sph_Ik_),
_precalcSH_(solvin->_precalcSH_),
mu_(solvin->mu_)
{
}


void Solver::precalc_sh_lf_lh()
{
  Pt q;
  vector<int> exp_pts;
  shared_ptr<BaseMolecule> mol;
  for (int i = 0; i < _sys_->get_n(); i++)
  {
    mol = _sys_->get_moli(i);
    for (int k = 0; k < mol->get_ns(); k++)
    {
      exp_pts = mol->get_gdpt_expj(k);
      for (int h=0; h < exp_pts.size(); h++)
      {
        q = mol->get_gridjh(k, exp_pts[h]);
        _precalcSH_->calc_and_add(q, _shCalc_);
      }
    }
  }
}

void Solver::precalc_sh_numeric()
{
  vector<int> exp_pts;
  Pt loc, sph_dist;
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
        for (int l = 0; l < _sys_->get_Ns_i(I); l++)
        {
          if (k == l) continue;
          else if (_T_->is_analytic(I, k, I, l)) continue;
          else
          {
            exp_pts = _sys_->get_gdpt_expij(I, l);
            for (int h = 0; h < exp_pts.size(); h++)
            {
              sph_dist = _sys_->get_centerik(I, k) - _sys_->get_centerik(I, l);
              loc = _sys_->get_gridijh(I, l, exp_pts[h]) - sph_dist;
              _precalcSH_->calc_and_add(loc, _shCalc_);
            }
          }
        }
    }
  }
}

void Solver::set_H_F(vector<shared_ptr<HMatrix > > h_spol,
                     vector<shared_ptr<FMatrix > > f_spol)
{
  int molt;
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    molt = _sys_->get_moli(I)->get_type();
    _H_[I]->set_all_mats(h_spol[molt]);
    _F_[I]->set_all_mats(f_spol[molt]);
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
      hnm_prev = ((inner) ? _prevH_[I]->get_mat_knm(k, n, m) :
                            _outerH_[I]->get_mat_knm(k, n, m));
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


void Solver::iter_innerH(int I, int k)
{
  double dev(1e20), toli(1e-6);
  int ct(0), mxCt(20);

  auto mol = _sys_->get_moli(I);
  while ( dev > toli && ct < mxCt )
  {
    _F_[I]->calc_vals(mol,_F_[I],_XF_[I],_H_[I],_IE_[I],_bCalc_,k,kappa_);
    _H_[I]->calc_vals(mol,_H_[I],_XH_[I],_F_[I],_IE_[I],_bCalc_,k);

    dev = calc_converge_H(I,k,true);
    update_prevH(I, k);

    ct++;
  }
}

void Solver::step(int t, int I, int k)
{
  // Do full step for spol and for mpol at t > 0
  if ((_sys_->get_n() == 1) || (t > 0))
  {
    _LH_[I]->init(_sys_->get_moli(I),_H_[I], _shCalc_, _bCalc_,
                  _precalcSH_, _expConsts_);
    _LH_[I]->calc_vals(_T_, _H_[I], _precalcSH_, k);

    _LF_[I]->init(_sys_->get_moli(I), _F_[I],_shCalc_, _bCalc_,
                  _precalcSH_, _expConsts_);
    _LF_[I]->calc_vals(_T_, _F_[I], _sys_, _precalcSH_, k);

    if (_sys_->get_n()>1)
      _LHN_[I]->calc_vals(_sys_, _T_, _rotH_, k); // use rotated H for mpol
    else
      _LHN_[I]->calc_vals(_sys_, _T_, _H_, k);
  }

  _XF_[I]->calc_vals(_sys_->get_moli(I), _bCalc_,
                     _LH_[I], _LF_[I], _LHN_[I], kappa_, k);
  _XH_[I]->calc_vals(_sys_->get_moli(I), _bCalc_, _LH_[I],
                     _LF_[I], _LHN_[I], kappa_, k);
}

double Solver::iter(int t)
{
  double mu_int(0), mu_mpol(0);
  if ( t == 0 ) mu_ = 0.0;
  Ns_tot_ = 0;

  if ((_sys_->get_n()>1) && (t==0)) // For 1st step of mpol
    for (int I = 0; I < _sys_->get_n(); I++)
    {
      auto molI = _sys_->get_moli(I);
      for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      {
        _LH_[I]->init(molI,_H_[I], _shCalc_, _bCalc_, _precalcSH_, _expConsts_);
        _LH_[I]->calc_vals(_T_, _H_[I], _precalcSH_, k);

        _LF_[I]->init(molI,_F_[I], _shCalc_, _bCalc_, _precalcSH_, _expConsts_);
        _LF_[I]->calc_vals(_T_, _F_[I], _sys_, _precalcSH_, k);

        _LHN_[I]->calc_vals(_sys_, _T_, _H_, k);
      //cout << "this is dev " << dev_sph_Ik_[I][k] << endl;
      }
    }

  for (int I = 0; I < _sys_->get_n(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      bool pol = true;
      // if there is more than 1 mol, run if there are sphs on mol to pol (LHN)
      if ((_sys_->get_n()>1) && (_LHN_[I]->get_interPol_k(k) != 0))
        pol = false;
      if((dev_sph_Ik_[I][k] > 0.1*mu_ && pol) ||
         ((t%5==0) && (_sys_->get_n()==1)))
      {
        update_outerH(I, k);
        step(t, I, k);
        iter_innerH(I, k);

        if ( t == 0 )
          update_rotH(I, k);

        dev_sph_Ik_[I][k] = calc_converge_H(I,k,false);
      //cout << "this is dev " << dev_sph_Ik_[I][k] << endl;
        mu_mpol += dev_sph_Ik_[I][k];
        Ns_tot_++;
        if ( dev_sph_Ik_[I][k] > mu_int )
          mu_int = dev_sph_Ik_[I][k];
      }
    } // end k
  }

  mu_ = mu_int;
  if (_sys_->get_n()==1)
    return mu_int;
  else
    return (mu_mpol / (double) Ns_tot_);
}


void Solver::update_prev_all()
{
  for (int I = 0; I < _H_.size(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      update_prevH(I, k);
      update_outerH(I, k);
    }
  }
}

void Solver::update_LHN_all()
{
  for (int I = 0; I < _H_.size(); I++)
  {
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      _LHN_[I]->calc_vals(_sys_, _T_, _H_, k);
    }
  }
}

void Solver::update_outerH(int I, int k)
{

  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
      _outerH_[I]->set_mat_knm(k, n, m, _H_[I]->get_mat_knm(k, n, m));
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
    }
  }
}

void Solver::update_rotH(int I, int k)
{

  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
      _rotH_[I]->set_mat_knm(k, n, m, _H_[I]->get_mat_knm(k, n, m));
    }
  }
}

void Solver::solve(double tol, int maxiter)
{
  double mu(1e15);
  for (int t = 0; t < maxiter; t++)
  {
    if (((t+1)%50) == 0) cout << "this is t " << t << endl;
    mu = iter(t);
    if (mu < tol) break;
  }

//for (int I = 0; I < _sys_->get_n(); I++)
//{
//  for (int k = 0; k < _sys_->get_Ns_i(I); k++)
//  {
//    cout << "This is LHN" << endl;
//    _LHN_[I]->print_kmat(k);
//  }
//}
//for (int I = 0; I < _sys_->get_n(); I++)
//{
//  for (int k = 0; k < _sys_->get_Ns_i(I); k++)
//  {
//    cout << "This is H" << endl;
//    _H_[I]->print_kmat(k);
//  }
//}
}


void Solver::reset_all()
{
  for (int I = 0; I < _sys_->get_n(); I++)
  {
    _E_[I]->calc_vals(_sys_->get_moli(I), _shCalc_,
                      _consts_->get_dielectric_prot());
    _LE_[I]->calc_vals(_sys_->get_moli(I), _shCalc_,
                       _consts_->get_dielectric_prot());
    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      _LF_[I]->reset_mat(k);
      _LH_[I]->reset_mat(k);
      _LHN_[I]->reset_mat(k);
      _XF_[I]->reset_mat(k);
      _XH_[I]->reset_mat(k);
      _prevH_[I]->reset_mat(k);
      _outerH_[I]->reset_mat(k);
      _rotH_[I]->reset_mat(k);
    }
  }

  _precalcSH_->clear_sh();
  precalc_sh_lf_lh();
  precalc_sh_numeric();
}


GradSolver::GradSolver(shared_ptr<SystemSAM> _sys,
                       shared_ptr<Constants> _consts,
                       shared_ptr<SHCalc> _shCalc,
                       shared_ptr<BesselCalc> _bCalc,
                       shared_ptr<TMatrix> _T,
                       vector<shared_ptr<FMatrix> > _F,
                       vector<shared_ptr<HMatrix> > _H,
                       vector<shared_ptr<IEMatrix> > _IE,
                       vector<vector<int > > interpol,
                       shared_ptr<PreCalcSH> precalc_sh,
                       shared_ptr<ExpansionConstants> _expConst,
                       int p, bool no_pre_sh)
:p_(p), Ns_tot_(0), _F_(_F), _H_(_H), _T_(_T),
_bCalc_(_bCalc), _shCalc_(_shCalc),
_sys_(_sys), _consts_(_consts), kappa_(_consts->get_kappa()),
interpol_(interpol), _IE_(_IE), _expConsts_(_expConst),
dev_sph_Ik_(_sys->get_n()),
dF_(_sys->get_n(), vector<shared_ptr<GradFMatrix> > (_sys->get_n())),
dH_(_sys->get_n(), vector<shared_ptr<GradHMatrix> > (_sys->get_n())),
prev_dH_(_sys->get_n(), vector<shared_ptr<GradHMatrix> > (_sys->get_n())),
outer_dH_(_sys->get_n(), vector<shared_ptr<GradHMatrix> > (_sys->get_n())),
dWF_(_sys->get_n(), vector<shared_ptr<GradWFMatrix> > (_sys->get_n())),
dWH_(_sys->get_n(), vector<shared_ptr<GradWHMatrix> > (_sys->get_n())),
dLF_(_sys->get_n(), vector<shared_ptr<GradLFMatrix> > (_sys->get_n())),
dLH_(_sys->get_n(), vector<shared_ptr<GradLHMatrix> > (_sys->get_n())),
dLHN_(_sys->get_n(), vector<shared_ptr<GradLHNMatrix> > (_sys->get_n())),
gradT_A_(_sys->get_n(), vector<shared_ptr<GradCmplxMolMat> > (_sys->get_n())),
precalcSH_(precalc_sh), noPreSH_(no_pre_sh)
{
  dF_.reserve(_sys_->get_n());
  for (int I = 0; I < _sys_->get_n(); I++) // With respect to
  {
    for (int J = 0; J < _sys_->get_n(); J++) // molecule
    {
      dF_[I][J] = make_shared<GradFMatrix> (J, I, _sys_->get_Ns_i(I), p_);
      dWF_[I][J] = make_shared<GradWFMatrix> (J, I, _sys_->get_Ns_i(I), p_,
                                              _consts->get_dielectric_prot(),
                                              _consts->get_dielectric_water(),
                                              _consts->get_kappa());

      dLF_[I][J] = make_shared<GradLFMatrix> (J, I, _sys_->get_Ns_i(I), p_);
      dLHN_[I][J] = make_shared<GradLHNMatrix> (J, I, _sys_->get_Ns_i(I), p_);

      prev_dH_[I][J] = make_shared<GradHMatrix> (J, I, _sys_->get_Ns_i(I),
                                                 p_, kappa_);
      outer_dH_[I][J] = make_shared<GradHMatrix> (J, I, _sys_->get_Ns_i(I),
                                                  p_, kappa_);
      dH_[I][J] = make_shared<GradHMatrix> (J,I,_sys_->get_Ns_i(I),p_,kappa_);

      dWH_[I][J] = make_shared<GradWHMatrix> (J,I,_sys_->get_Ns_i(I),p_,kappa_);
      dLH_[I][J] = make_shared<GradLHMatrix> (J,I,_sys_->get_Ns_i(I),p_,kappa_);

      gradT_A_[I][J] = make_shared<GradCmplxMolMat> (J,I,_sys_->get_Ns_i(I),
                                                  p_);
    }
    dev_sph_Ik_[I].resize(_sys_->get_Ns_i(I));
  }

  for (int i = 0; i < _T_->get_T_ct(); i++)  _T_->compute_derivatives_i(i);
}

GradSolver::GradSolver(shared_ptr<GradSolver> gradin)
: p_(gradin->p_), kappa_(gradin->kappa_), Ns_tot_(gradin->Ns_tot_),
precalcSH_(gradin->precalcSH_),
noPreSH_(gradin->noPreSH_),
_F_(gradin->_F_),
_H_(gradin->_H_),
_IE_(gradin->_IE_),
_T_(gradin->_T_),
dF_(gradin->dF_),
dH_(gradin->dH_),
prev_dH_(gradin->prev_dH_),
outer_dH_(gradin->outer_dH_),
dWF_(gradin->dWF_),
dWH_(gradin->dWH_),
dLF_(gradin->dLF_),
dLH_(gradin->dLH_),
dLHN_(gradin->dLHN_),
gradT_A_(gradin->gradT_A_),
_sys_(gradin->_sys_),
_shCalc_(gradin->_shCalc_),
_bCalc_(gradin->_bCalc_),
_expConsts_(gradin->_expConsts_),
_consts_(gradin->_consts_),
interpol_(gradin->interpol_),
dev_sph_Ik_(gradin->dev_sph_Ik_)
{ }

void GradSolver::solve(double tol, int maxiter)
{
  double mu;
  double scale_dev = (double)(p_*(p_+1)*0.5*3.0);
  int j, ct;

  pre_compute_gradT_A();

  for ( j = 0; j < _sys_->get_n(); j++ ) // gradient WRT j
  {
    ct = 0;
    mu = scale_dev;
    while(mu > tol)
    {
      mu = iter(ct, j);
      if ((ct+1) % 10 == 0) cout << "Iter step " << ct << endl;
      if (ct > maxiter)  break;
      ct++;
    }
  }
}

double GradSolver::iter(int t, int wrt)
{
  double inter_pol_d(10.), mu_mpol(0);
  Ns_tot_ = 0;

  shared_ptr<BaseMolecule> molI;
  vector<double> besseli, besselk;
  vector<int> mol_loop;

  for (int I = 0; I < _sys_->get_n(); I++)
    if (( I != wrt ) && ( _sys_->get_min_dist(wrt, I) < inter_pol_d))
      mol_loop.push_back(I);

  mol_loop.push_back(wrt);

  for (int ind = 0; ind < mol_loop.size(); ind++)
  {
    int I = mol_loop[ind];
    molI = _sys_->get_moli(I);
    if ((_sys_->get_min_dist(wrt, I) > inter_pol_d) && (I != wrt)) continue;

    dLHN_[wrt][I]->calc_all_vals(_sys_, _T_, gradT_A_[wrt], dH_[wrt]);

    for (int k = 0; k < _sys_->get_Ns_i(I); k++)
    {
      if (interpol_[I][k] != 0) continue;
      if (!_sys_->get_moli(I)->is_J_in_interk(k, wrt) && (I != wrt))
        continue;

      besseli = _bCalc_->calc_mbfI(p_+1, kappa_*molI->get_ak(k));
      besselk = _bCalc_->calc_mbfK(p_+1, kappa_*molI->get_ak(k));

      update_outer_gradH(I, wrt, k);
      step( t, I, wrt, k, besseli, besselk);
      iter_inner_gradH(I, wrt, k, besseli, besselk);

      dLF_[wrt][I]->init_k(k,molI,dF_[wrt][I],_shCalc_, precalcSH_ ,
                           _expConsts_, noPreSH_);
      dLH_[wrt][I]->init_k(k,molI,dH_[wrt][I],_shCalc_,_bCalc_, precalcSH_,
                           _expConsts_, noPreSH_);

      dev_sph_Ik_[I][k] = calc_converge_gradH(I, wrt, k, false);
      mu_mpol += dev_sph_Ik_[I][k];
      Ns_tot_++;
    }
  }
  return (mu_mpol / (double) Ns_tot_);
}

// Perform updates of various expansions
void GradSolver::step(int t, int I, int wrt, int k, vector<double> &besseli,
                      vector<double> &besselk)
{
  shared_ptr<BaseMolecule> molI = _sys_->get_moli(I);
  dLF_[wrt][I]->calc_val_k(k, molI, interpol_[I], _T_, dF_[wrt][I],
                           precalcSH_, noPreSH_);
  dLH_[wrt][I]->calc_val_k(k, molI, interpol_[I], _T_, dH_[wrt][I],
                           precalcSH_, noPreSH_);
}


void GradSolver::iter_inner_gradH(int I, int wrt, int k,
                                  vector<double> &besseli,
                                  vector<double> &besselk)
{
  double dev(1e20), toli(1e-6);
  int ct(0), mxCt(20);

  auto mol = _sys_->get_moli(I);
  while ( dev > toli && ct < mxCt )
  {
    dWF_[wrt][I]->calc_val_k(k,mol,besseli,besselk,dH_[wrt][I],dF_[wrt][I],
                             dLH_[wrt][I], dLHN_[wrt][I], dLF_[wrt][I]);
    dF_[wrt][I]->calc_val_k(k, _IE_[I], dWF_[wrt][I]);
    dWH_[wrt][I]->calc_val_k(k, mol, besseli, besselk, dH_[wrt][I],
                             dF_[wrt][I], dLH_[wrt][I], dLHN_[wrt][I],
                             dLF_[wrt][I]);
    dH_[wrt][I]->calc_val_k(k, besseli, _IE_[I], dWH_[wrt][I]);

    dev = calc_converge_gradH(I, wrt, k, true);
    update_prev_gradH(I, wrt, k);

    ct++;
  }
}


double GradSolver::calc_converge_gradH(int I, int wrt, int k, bool inner)
{
  double mu(0), rt(0), num(0), den(0);
  cmplx hnm_curr, hnm_prev;
  double PRC_LM = 1e-30;

  for (int d = 0; d < 3; d++)
  {
    for (int n = 0; n < dH_[wrt][I]->get_p(); n++)
    {
      for (int m = -n; m < n+1; m++)
      {
        hnm_curr = dH_[wrt][I]->get_mat_knm_d(k, n, m, d);
        hnm_prev = ((inner) ? prev_dH_[wrt][I]->get_mat_knm_d(k, n, m, d) :
                    outer_dH_[wrt][I]->get_mat_knm_d(k, n, m, d));
        if (inner)
        {
          if ( m < 0 ) continue;
          if ((fabs(hnm_prev.real())<PRC_LM)&&(fabs(hnm_prev.imag())<PRC_LM))
          {
            if ((fabs(hnm_curr.real())<PRC_LM)&&(fabs(hnm_curr.imag())<PRC_LM))
              continue;
            else
              rt = 1.0;
          } else
          {
            if ((fabs(hnm_curr.real())<PRC_LM)&&(fabs(hnm_curr.imag())<PRC_LM))
              rt = 1.0;
            else
              rt = 2.0*norm(hnm_curr-hnm_prev)/(norm(hnm_curr)+norm(hnm_prev));
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
  }

  if (inner)
    return mu / double(4.0*p_*p_);
  else
  {
    if ( fabs(den) > 1e-15 )
      mu = (num / den);
    else
      mu = 0.0;
    return mu / 4.0 / 3.0;
  }
}

void GradSolver::update_outer_gradH(int I, int wrt, int k)
{

  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
      outer_dH_[wrt][I]->set_mat_knm(k, n, m,
                                     dH_[wrt][I]->get_mat_knm(k, n, m));
    }
  }
}

void GradSolver::update_prev_gradH(int I, int wrt, int k)
{

  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m <= n; m++)
    {
      prev_dH_[wrt][I]->set_mat_knm(k, n, m, dH_[wrt][I]->get_mat_knm(k, n, m));
    }
  }
}

// precompute gradT times H(i,j) for all pairs of MoleculeSAMs
void GradSolver::pre_compute_gradT_A()
{
  shared_ptr<BaseMolecule> molI;
  double aIk, aJl, dist, cut_act(100.0), cut_er(10.0);
  Pt Ik, Jl, v;
  MyMatrix<Ptx> reex(p_, 2*p_+1);

  for (int I = 0; I < _sys_->get_n(); I++)
    for (int J = 0; J < _sys_->get_n(); J++)
      gradT_A_[J][I]->reset_mat();

  for (int J = 0; J < _sys_->get_n(); J++)
  {
  for (int I = 0; I < _sys_->get_n(); I++)
  {
      if ( I == J ) continue;
      for (int l = 0; l < _sys_->get_Ns_i(J); l++)
      {
        Jl = _sys_->get_centerik(J, l);
        aJl = _sys_->get_aik(J, l);
      for (int k = 0; k < _sys_->get_Ns_i(I); k++)
      {
        Ik = _sys_->get_centerik(I, k);
        aIk = _sys_->get_aik(I, k);

          dist = _sys_->get_pbc_dist_vec_base(Ik, Jl).norm();
          if ( dist < (cut_er+aIk+aJl) )
          {
            reex = _T_->re_expandX_gradT(_H_[J]->get_mat_k(l), I, k, J, l);
            gradT_A_[J][I]->add_mat_k(k, reex);

            if (interpol_[J][l] == 0)
            {
              reex = _T_->re_expandX_gradT(_H_[I]->get_mat_k(k), J, l, I, k);
              gradT_A_[J][J]->add_mat_k(l, reex);
            }
          } else if ((dist < (cut_act+aIk+aJl)) && (interpol_[J][l] == 0))
          {

            reex = _T_->re_expandX_gradT(_H_[I]->get_mat_k(k), J, l, I, k);
            gradT_A_[J][J]->add_mat_k(l, reex);


          } else if ((interpol_[J][l] != 0) &&
                     (_sys_->get_moli(J)->get_inter_act_k(l).size() != 0) &&
                     (dist < (cut_act+aIk+aJl)))
          {
//            cout << "Reex 100A for mol (org) " << I << " sph " << k
//            << " to dest : " << J << " and sph " << l  << endl;
            reex = _T_->re_expandX_gradT(_H_[I]->get_mat_k(k), J, l, I, k);
            gradT_A_[J][J]->add_mat_k(l, reex);
          }
        } // end l
      }
    }
  }
}

} /* namespace pbsolvers */

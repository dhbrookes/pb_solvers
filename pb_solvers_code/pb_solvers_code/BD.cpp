//
//  BD.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BD.h"
#include <iostream>


BDStep::BDStep(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
       vector<double> trans_diff_consts,
       vector<double> rot_diff_consts,
       bool diff, bool force)
:transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

BDStep::BDStep(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
               bool diff, bool force)
:transDiffConsts_(_sys->get_n()), rotDiffConsts_(_sys->get_n()),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  for (int i = 0; i < _sys_->get_n(); i++)
  {
    transDiffConsts_[i] = _sys_->get_dtransi(i);
    rotDiffConsts_[i] = _sys_->get_droti(i);
  }
  
  random_device rd;
  randGen_ = mt19937(rd());
}


double BDStep::compute_dt( double dist )
{
  double DISTCUTOFF_TIME = 20.0;
  if ( dist - DISTCUTOFF_TIME > 0 )
    return 2.0 + ( dist - DISTCUTOFF_TIME )/15.0;
  else
    return 2.0;
}

double BDStep::compute_min_dist( )
{
  int i, j;
  Pt pt1, pt2;
  double newD, minDist = 1e100;
  for (i = 0; i < _sys_->get_n(); i++)
  {
    for (j = i+1; j < _sys_->get_n(); j++)
    {
      pt1 = _sys_->get_centeri(i); pt2 = _sys_->get_centeri(j);
      newD = pt1.dist(pt2);
      if (newD < minDist) minDist = newD;
    }
  }
  return minDist;
}

Pt BDStep::rand_vec(double mean, double var)
{
  normal_distribution<double> dist (mean, sqrt(var));
  Pt pout = Pt(dist(randGen_), dist(randGen_), dist(randGen_));
  return pout;
}


bool BDStep::check_for_collision(int mol, Pt new_pt)
{
  bool collision = false;
  int j;
  double dist, aj;
  Pt pj;
  double ai = _sys_->get_ai(mol);
  
  for (j = 0; j < _sys_->get_n(); j++)
  {
    if (j == mol) continue;
    pj = _sys_->get_centeri(j);
    aj = _sys_->get_ai(j);
    
    dist = new_pt.dist(pj);
    if (dist < (ai + aj))
    {
      collision = true;
      break;
    }
  }
  return collision;
}

void BDStep::indi_trans_update(int i, MyVector<double> fi)
{
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = transDiffConsts_[i] * dt_ * ikT_int;
  Pt dr = Pt(fi * coeff);
  Pt center = _sys_->get_centeri(i);
  Pt rand, new_pt;
  bool accept = false;
  
  int tries = 0;
  while (!accept && tries < 500)
  {
    tries++;
    rand = (diff_) ? rand_vec(0,2*transDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    new_pt = center + (dr + rand);
    if (!check_for_collision(i, new_pt))
    {
      _sys_->translate_mol(i, dr+rand);
      accept = true;
    }
  }
}


void BDStep::indi_rot_update(int i, MyVector<double> tau_i)
{
  Pt rand, new_pt;
  Quat qrot;
  
  bool accept = false;
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = rotDiffConsts_[i] * dt_ * ikT_int;
  
  Pt dtheta = (tau_i * coeff);
  Pt center = _sys_->get_centeri(i);
  
  while (! accept)
  {
    rand = (diff_) ? rand_vec(0,2*rotDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    dtheta = dtheta + rand;
    
    // creating zero quaternion if there is no rot
    if (dtheta.norm() < 1e-15) qrot = Quat();
    else qrot = Quat(dtheta.norm(), dtheta);
    new_pt = qrot.rotate_point(center);
    
    if (!check_for_collision(i, new_pt))
    {
      _sys_->rotate_mol(i, qrot);
      accept = true;
    }
  }
}

void BDStep::bd_update(shared_ptr<VecOfVecs<double>::type> _F,
                   shared_ptr<VecOfVecs<double>::type> _tau)
{
  int i;
  dt_ = compute_dt(compute_min_dist());
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( transDiffConsts_[i] != 0) indi_trans_update(i, _F->operator[](i));
  }
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( rotDiffConsts_[i] != 0) indi_rot_update(i, _tau->operator[](i));
  }
  update_sys_time(dt_);
}


BDRun::BDRun(shared_ptr<ASolver> _asolv,
             shared_ptr<BaseTerminate> _terminator, shared_ptr<System> _sys,
             shared_ptr<Constants> _consts,
             bool diff, bool force, int maxiter, double prec)
:maxIter_(maxiter), _asolver_(_asolv), prec_(prec), _terminator_(_terminator)
{
  _fCalc_ = make_shared<ForceCalc>(_asolv);
  _torCalc_ = make_shared<TorqueCalc>(_asolv);
  _stepper_ = make_shared<BDStep> (_sys, _consts, diff, force);
  
}


void BDRun::run()
{
  
  int i = 0;
  bool term = false;
  while (i < maxIter_ and !term)
  {
    _asolver_->reset_all(_stepper_->get_system());
    _asolver_->solve_A(prec_);
    _asolver_->solve_gradA(prec_);
    
    _fCalc_->calc_force();
    _torCalc_->calc_tau();
    
    _stepper_->bd_update(_fCalc_->get_F(), _torCalc_->get_Tau());
    if (_terminator_->is_terminated(_stepper_->get_system()))
    {
      term = true;
    }
    i++;
  }
  
  
}
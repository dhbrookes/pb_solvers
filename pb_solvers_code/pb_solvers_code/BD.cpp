//
//  BD.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BD.h"
#include <iostream>


BD::BD(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
       vector<double> trans_diff_consts,
       vector<double> rot_diff_consts,
       bool diff, bool force)
:transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

double BD::compute_dt( double dist )
{
  double DISTCUTOFF_TIME = 20.0;
  if ( dist - DISTCUTOFF_TIME > 0 )
    return 2.0 + ( dist - DISTCUTOFF_TIME )/15.0;
  else
    return 2.0;
}

double BD::compute_min_dist( )
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

Pt BD::rand_vec(double mean, double var)
{
  normal_distribution<double> dist (mean, sqrt(var));
  Pt pout = Pt(dist(randGen_), dist(randGen_), dist(randGen_));
  return pout;
}


bool BD::check_for_collision(int mol, Pt new_pt)
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

void BD::indi_trans_update(int i, MyVector<double> fi)
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


void BD::indi_rot_update(int i, MyVector<double> tau_i)
{
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = transDiffConsts_[i] * dt_ * ikT_int;
  Pt dtheta = (tau_i * coeff);
  double f0 = (abs(tau_i[0])<1e-15) ? 0:tau_i[0];
  double f1 = (abs(tau_i[1])<1e-15) ? 0:tau_i[1];
  double f2 = (abs(tau_i[2])<1e-15) ? 0:tau_i[2];
  cout << " Mol " << i << " Dtr " << rotDiffConsts_[i] <<" dt " << dt_ << " & ikt " << ikT_int <<  " Taui " << f0 << "  " << f1 << "  " << f2 << endl;
  bool accept = false;
  Pt center = _sys_->get_centeri(i);
  Pt rand, new_pt;
  while (! accept)
  {
    rand = (diff_) ? rand_vec(0,2*rotDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    dtheta = dtheta + rand;
    Quat qrot (dtheta.norm(), dtheta.x(), dtheta.y(), dtheta.z());
    cout << " This is my quat " << qrot.get_w() << " a " << qrot.get_a() << " b " << qrot.get_b() << " c " << qrot.get_c() << endl;
    
    new_pt = qrot.rotate_point(center);
    if (!check_for_collision(i, new_pt))
    {
      _sys_->rotate_mol(i, qrot);
      accept = true;
    }
  }
}

void BD::bd_update(VecOfVecs<double>::type F, VecOfVecs<double>::type tau)
{
  int i;
  dt_ = compute_dt(compute_min_dist());
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( transDiffConsts_[i] != 0) indi_trans_update(i, F[i]);
  }
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( rotDiffConsts_[i] != 0) indi_rot_update(i, tau[i]);
  }
  update_sys_time(dt_);
}

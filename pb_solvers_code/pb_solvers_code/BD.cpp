//
//  BD.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BD.h"


BD::BD(System sys, double dt, vector<double> trans_diff_consts,
   vector<double> rot_diff_consts)
:dt_(dt), transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts)
{
  _sys_ = make_shared<System>(sys);
  
  random_device rd;
  randGen_ = mt19937(rd());
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
  double coeff = transDiffConsts_[i] * dt_ * _sys_->get_consts().get_ikbt();
  Pt dr = Pt(fi * coeff);
  bool accept = false;
  Pt center = _sys_->get_centeri(i);
  Pt rand, new_pt;
  while (!accept)
  {
    rand = rand_vec(0, 2 * rotDiffConsts_[i] *dt_);
    new_pt = center + (dr + rand);
    if (!check_for_collision(i, new_pt))
    {
      _sys_->translate_mol(i, dr);
      accept = true;
    }
  }
}


void BD::indi_rot_update(int i, MyVector<double> tau_i)
{
  double coeff = (rotDiffConsts_[i] * dt_) / (kb_ * T_);
  Pt dtheta = (tau_i * coeff);
  
  bool accept = false;
  Pt center = _sys_->get_centeri(i);
  Pt rand, new_pt;
  while (! accept)
  {
    rand = rand_vec(0, 2 * transDiffConsts_[i] * dt_);
    dtheta = dtheta + rand;
    Quat qrot (dtheta.norm(), dtheta.x(), dtheta.y(), dtheta.z());
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
  for (i = 0; i < _sys_->get_n(); i++)
  {
    indi_trans_update(i, F[i]);
  }
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    indi_rot_update(i, tau[i]);
  }
  
  
}

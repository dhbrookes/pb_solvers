//
//  BD.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BD.h"


BD::BD(System sys, double dt, vector<double> trans_diff_consts,
   vector<double> rot_diff_consts, double kb, double T)
:dt_(dt), transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts),
kb_(kb), T_(T)
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


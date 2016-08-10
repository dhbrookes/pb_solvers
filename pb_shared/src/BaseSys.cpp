//
//  BaseSys.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 8/9/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BaseSys.h"



BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, double drot,
                           double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(type_idx),
Nc_((int) pos.size()), centers_(1), as_(1), Ns_(1)
{
  
  set_Dtr_Drot(moveType_);
}

BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, vector<Pt> cens,
                           vector<double> as, double drot, double dtrans)
:BaseMolecule(type, type_idx, movetype, qs, pos, vdwr, drot, dtrans)
{
  centers_ = cens;
  as_ = as;
  Ns_ = (int) centers_.size();
  
  
  set_Dtr_Drot(moveType_);
}

BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, Pt cen,
                           double a, double drot, double dtrans)
:BaseMolecule(type, type_idx, movetype, qs, pos, vdwr, drot, dtrans)
{
  centers_[0] = cen;
  as_[0] = a;
  set_Dtr_Drot(moveType_);
}

void BaseMolecule::set_Dtr_Drot(string type)
{
  if ((type == "stat") or (type == "rot"))  dtrans_ = 0.0;
  if (type == "stat") drot_ = 0.0;
}

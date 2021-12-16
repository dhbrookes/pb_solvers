//
//  BaseBD.cpp
//  pbam_xcode
//
//  Created by David Brookes on 9/1/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BaseBD.h"

namespace pbsolvers
{

string BaseTerminate::get_how_term(shared_ptr<BaseSystem> _sys)
{
  return "";
}

const bool BaseTerminate::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  return false;
}


TimeTerminate::TimeTerminate(double end_time)
:BaseTerminate(), endTime_(end_time)
{
  char buff[400];
  sprintf(buff, "System has fulfilled condition: time >= %7.1f;\t", endTime_);
  how_term_ = buff;
}

const bool TimeTerminate::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  bool done = false;
  if (_sys->get_time() >= endTime_) done = true;
  return done;
}


CoordTerminate::CoordTerminate(int mol_idx, CoordType coord_type,
                               BoundaryType bound_type, double bound_val)
:BaseTerminate(), molIdx_(mol_idx), coordType_(coord_type),
boundType_(bound_type), boundaryVal_(bound_val)
{
  char buff[400];
  string cord = "r", eq = ">=";
  if (coordType_ == X)      cord = "x";
  else if (coordType_ == Y) cord = "y";
  else if (coordType_ == Z) cord = "z";
  
  if (boundType_ == LEQ)    eq   = "<=";
  
  sprintf(buff, "MoleculeAM type %d has fulfilled condition: %s %s %5.2f;\t",
          molIdx_, cord.c_str(), eq.c_str(), boundaryVal_);
  how_term_ = buff;
}


const bool CoordTerminate::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  bool done = false;
  int i, idx;
  for ( i = 0; i < _sys->get_typect(molIdx_); i++)
  {
    idx = _sys->get_mol_global_idx( molIdx_, i);
    Pt mol_coord = _sys->get_unwrapped_center(idx);
    double test_val;
    if (coordType_ == X)      test_val = mol_coord.x();
    else if (coordType_ == Y) test_val = mol_coord.y();
    else if (coordType_ == Z) test_val = mol_coord.z();
    else                      test_val = mol_coord.norm();
    
    if ((boundType_ == LEQ) && (test_val <= boundaryVal_))      return true;
    else if ((boundType_ == GEQ) && (test_val >= boundaryVal_)) return true;
  }
  return done;
}


CombineTerminate::CombineTerminate(vector<shared_ptr<BaseTerminate> > terms,
                                   HowTermCombine how_combine)
:BaseTerminate(), terms_(terms), howCombine_(how_combine)
{
}

const bool CombineTerminate::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  bool done;
  howCombine_== ALL ? done=true : done=false;
  for (int i = 0; i < terms_.size(); i++)
  {
    if (terms_[i]->is_terminated(_sys) == ! done)
    {
      done=!done;
      break;
    }
  }
  return done;
}

string CombineTerminate::get_how_term(shared_ptr<BaseSystem> _sys)
{
  string how_term = "";
  bool done;
  howCombine_== ALL ? done=true : done=false;
  for (int i = 0; i < terms_.size(); i++)
  {
    if (terms_[i]->is_terminated(_sys) == true)
    {
      how_term += terms_[i]->get_how_term(_sys);
    }
  }
  return how_term;
}

BaseBDStep::BaseBDStep(shared_ptr<BaseSystem> _sys,
                       shared_ptr<Constants> _consts,
                       vector<double> trans_diff_consts,
                       vector<double> rot_diff_consts,
                       bool diff, bool force)
:transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

BaseBDStep::BaseBDStep(shared_ptr<BaseSystem> _sys,
                       shared_ptr<Constants> _consts,
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


double BaseBDStep::compute_dt( )
{
  double DISTCUTOFF_TIME = 20.0;
  if ( min_dist_ - DISTCUTOFF_TIME > 0 )
  return 2.0 + ( min_dist_ - DISTCUTOFF_TIME )/15.0;
  else
  return 2.0;
}

Pt BaseBDStep::rand_vec(double mean, double var)
{
  normal_distribution<double> dist (mean, sqrt(var));
  Pt pout = Pt(dist(randGen_), dist(randGen_), dist(randGen_));
  return pout;
}

void BaseBDStep::indi_trans_update(int i, Pt fi)
{
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = transDiffConsts_[i] * dt_ * ikT_int;

  Pt dr = Pt(fi * coeff);
  Pt rand, new_pt;
  bool accept = false;
  
  int tries = 0;
  while (!accept && tries < 500)
  {
    tries++;
    rand = (diff_) ? rand_vec(0,2*transDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    _sys_->translate_mol(i, dr + rand);
    accept = true;
    try {
      _sys_->check_for_overlap();
    } catch (OverlappingMoleculeException) {
      _sys_->translate_mol(i, (dr+rand) * -1);
      accept = false;
    }
  }
}


void BaseBDStep::indi_rot_update(int i, Pt tau_i)
{
  Pt rand, new_pt;
  Quat qrot;
  
  bool accept = false;
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = rotDiffConsts_[i] * dt_ * ikT_int;
  
  Pt dtheta = (tau_i * coeff);
  
  while (! accept)
  {
    rand = (diff_) ? rand_vec(0,2*rotDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    dtheta = dtheta + rand;
    
    // creating zero quaternion if there is no rot
    if (dtheta.norm() < 1e-15) qrot = Quat();
    else qrot = Quat(dtheta.norm(), dtheta);
    _sys_->rotate_mol(i, qrot);
    accept = true;
    try {
      _sys_->check_for_overlap();
    } catch (OverlappingMoleculeException) {
      _sys_->rotate_mol(i, qrot.conj());
      accept = false;
    }
  }
}

void BaseBDStep::bd_update(shared_ptr<vector<Pt> > _F,
                           shared_ptr<vector<Pt> > _tau)
{
  int i;
  compute_min_dist();
  dt_ = compute_dt();
  
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

BaseBDRun::BaseBDRun(shared_ptr<BaseTerminate> _terminator, string outfname,
                     int num, bool diff, bool force, int maxiter, double prec)
:maxIter_(maxiter), prec_(prec), _terminator_(_terminator)
{
}

} /* namespace pbsolvers */

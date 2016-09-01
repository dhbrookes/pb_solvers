//
//  BaseBD.cpp
//  pbam_xcode
//
//  Created by David Brookes on 9/1/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BaseBD.h"


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




//
//  BaseBD.hpp
//  pbam_xcode
//
//  Created by David Brookes on 9/1/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BaseBD_h
#define BaseBD_h

#include <stdio.h>
#include <random>
#include <memory>
#include <vector>
#include "BaseSys.h"

using namespace std;
/*
 Base class for implementing termination conditions in BD
 */
class BaseTerminate
{
public:
  BaseTerminate() { }
  
  virtual const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  virtual string get_how_term(shared_ptr<BaseSystem> _sys);
};

enum CoordType { X, Y, Z, R };
enum BoundaryType { LEQ, GEQ };

/*
 Class for time based termination
 */
class TimeTerminate : public BaseTerminate
{
protected:
  double endTime_; //termination time
  string how_term_;
  
public:
  TimeTerminate(double end_time);
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
};

/*
 Class for coordinate based termination. This terminates based on whether
 the specified MoleculeAM satisfies the BoundaryType condition on the CoordType
 with the given boundary value.
 */
class CoordTerminate : public BaseTerminate
{
protected:
  double boundaryVal_;
  int molIdx_;
  CoordType coordType_;
  BoundaryType boundType_;
  string how_term_;
  
public:
  CoordTerminate(int mol_idx, CoordType coord_type,
                 BoundaryType bound_type, double bound_val);
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  
  string get_how_term(shared_ptr<BaseSystem> _sys)   { return how_term_; }
};

/*
 Combine termination conditions
 */
enum HowTermCombine { ALL, ONE };

class CombineTerminate: public BaseTerminate
{
protected:
  vector<shared_ptr<BaseTerminate> > terms_;
  HowTermCombine howCombine_;
  
public:
  CombineTerminate(vector<shared_ptr<BaseTerminate> > terms,
                   HowTermCombine how_combine);
  
  const bool is_terminated(shared_ptr<BaseSystem> _sys) const;
  
  string get_how_term(shared_ptr<BaseSystem> _sys);
};



#endif /* BaseBD_h */

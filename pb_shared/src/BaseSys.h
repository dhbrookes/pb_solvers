//
//  BaseSys.hpp
//  pbsam_xcode
//
//  Created by David Brookes on 8/9/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef BaseSys_h
#define BaseSys_h

#include <stdio.h>
#include <vector>
#include "util.h"

using namespace std;
/*
 Base class for storing relevant data about each MoleculeAM
 */
class BaseMolecule
{
protected:
  string              moveType_;
  int                 type_; // int index of type of MoleculeAM, 0 based
  int                 typeIdx_; // int index of mol within given type_, 0 based
  double              drot_;  // rotational diffusion coefficient
  double              dtrans_; // translational diffusion coefficients
  int                 Nc_;  // number of charges in this MoleculeAM
  
  vector<double>      qs_;  // magnitude of each charge in the MoleculeAM
  vector<Pt>          pos_;  // position of each charge in the MoleculeAM
  vector<double>      vdwr_; // van der waal radius of each charge
  
  int                 Ns_;  // number of coarse grained spheres
  vector<Pt>          centers_; //coarse-grained sphere centers
  vector<double>      as_; // coarse-grained sphere radii
  
  void set_Dtr_Drot(string type);
  
public:
  
  BaseMolecule()  { }
//  BaseMoleculeAM(const BaseMoleculeAM& mol);
  
  // no centers or raddi of cg spheres
  BaseMolecule(int type, int type_idx, string movetype, vector<double> qs,
               vector<Pt> pos, vector<double> vdwr, double drot=0,
               double dtrans=0);
  
  // list of cg centers and radii
  BaseMolecule(int type, int type_idx, string movetype, vector<double> qs,
               vector<Pt> pos, vector<double> vdwr, vector<Pt> cens,
               vector<double> as, double drot=0, double dtrans=0);
  
  // single center and radius of cg sphere
  BaseMolecule(int type, int type_idx, string movetype, vector<double> qs,
               vector<Pt> pos, vector<double> vdwr, Pt cen,
               double a, double drot=0, double dtrans=0);
  
  const int get_m() const               { return Nc_; } // called M in analytic
  const int get_nc() const              { return Nc_; }
  Pt get_posj(int j) const              { return pos_[j]; }
  string get_move_type() const          { return moveType_; }
  int get_type() const                  { return type_; }
  int get_type_idx() const              { return typeIdx_; }
  double get_drot() const               { return drot_; }
  double get_dtrans() const             { return dtrans_; }
  const double get_qj(int j) const      { return qs_[j]; }
  const double get_radj(int j) const    { return vdwr_[j]; }
  Pt get_posj_realspace(int j)          { return pos_[j] + get_cen_j(j); }
  const int get_ns() const              { return Ns_; }
  
  // get center for charge j
  virtual Pt get_cen_j(int j) = 0;
  
  virtual void translate(Pt dr, double boxlen) = 0;
  virtual void rotate(Quat qrot) = 0;
  virtual void rotate(MyMatrix<double> rotmat) = 0;
  
};


//class BaseSystem
//{
//  vector<BaseMolecule> mols;
//};

#endif /* BaseSys_h */

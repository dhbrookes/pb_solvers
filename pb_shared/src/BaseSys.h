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
#include <map>
#include <memory>
#include "util.h"
#include "Constants.h"

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
  const double get_ak(int k) const      { return as_[k]; }
  Pt get_centerk(int k) const           { return centers_[k]; }
  const double get_qj(int j) const      { return qs_[j]; }
  const double get_radj(int j) const    { return vdwr_[j]; }
  Pt get_posj_realspace(int j)          { return pos_[j] + get_cen_j(j); }
  const int get_ns() const              { return Ns_; }
  
  // get center for charge j
  virtual Pt get_cen_j(int j) = 0;
  
  virtual void translate(Pt dr, double boxlen) = 0;
  virtual void rotate(Quat qrot) = 0;
  virtual void rotate(MyMatrix<double> rotmat) = 0;
  
  
  //below are methods only for PBAM that should be overridden:
  virtual vector<int> get_pol()                 { return vector<int> (); }
  virtual vector<int> get_act()                 { return vector<int> (); }
  virtual void clear_inter_pol()                { }
  virtual void clear_inter_act()                { }
  virtual void add_J_to_pol(int J)              { }
  virtual void add_J_to_interact(int J)         { }
  virtual bool is_J_in_pol( int J )             { return true; }
  virtual bool is_J_in_interact( int J )        { return true; }
  
  
  //below are methods only for PBSAM that should be overridden
  virtual void set_gridj(int j, vector<Pt> grid) { }
  virtual void set_gridexpj(int j, vector<int> grid_exp) { }
  virtual void set_gridburj(int j, vector<int> grid_bur) { }
  virtual void add_Jl_to_interk( int k, int J, int l) { }
  virtual void add_Jl_to_inter_act_k( int k, int J, int l) { }
  virtual Pt get_gridjh(int j, int h) const   { return Pt(); }
  virtual Pt get_cog() const                  { return Pt(); }
  virtual Pt get_unwrapped_center() const     { return Pt(); }
  virtual vector<int> get_gdpt_expj(int j) const { return vector<int>(); }
  virtual int get_cg_of_ch(int k) const       { return 0; }
  virtual int get_nc_k(int k) const           { return 0; }
  virtual vector<int> get_neighj(int j) const { return vector<int> (); }
  virtual vector<Pt> get_gridj(int j) const   { return vector<Pt> (); }
  virtual vector<int> get_gdpt_burj(int j) const { return vector<int> (); }
  virtual vector<int> get_ch_allin_k(int k)   { return vector<int> (); }
  virtual vector<int> get_ch_allout_k(int k)  { return vector<int> (); }
  virtual int get_ch_k_alpha(int k, int alpha){ return 0; }
  virtual vector<vector<int> > get_inter_act_k(int k) {return vector<vector<int> >(); }
  virtual bool is_J_in_interk( int k, int J ) { return true; }
};


class BaseSystem
{
protected:
  
  int                                      N_;
  double                                   lambda_;
  vector<shared_ptr<BaseMolecule> >        molecules_;
  
  double                       boxLength_;
  double                       cutoff_;
  
  double                       t_;  // time in a BD simulatio
  int                          ntype_;  //count of unique molecules
  vector<int>                  typect_; //count of molecules of each
  
  map<vector<int>, int>        typeIdxToIdx_;
  
  const double calc_average_radius() const;
  
public:
  BaseSystem()  { }
  
  BaseSystem(vector<shared_ptr<BaseMolecule> > mols,
             double cutoff=Constants::FORCE_CUTOFF,
             double boxlength=Constants::MAX_DIST);
  
  const int get_n() const                  {return N_;}
  const int get_ntype()                    {return ntype_;}
  const int get_typect(int i)              {return typect_[i];}
  const double get_boxlength() const       {return boxLength_;}
  const double get_cutoff() const          {return cutoff_;}
  const double get_time() const            {return t_;}
  const double get_lambda() const          {return lambda_;}
  
  const int get_Nc_i(int i) const          { return molecules_[i]->get_nc();}
  const int get_Ns_i(int i) const          { return molecules_[i]->get_ns();}
  const double get_aik(int i, int k) const { return molecules_[i]->get_ak(k);}
  Pt get_centerik(int i, int k) const      { return molecules_[i]->get_centerk(k); }
  const double get_Mi(int i) const         {return molecules_[i]->get_m();}
  const double get_qij(int i, int j) const {return molecules_[i]->get_qj(j);}
  const double get_radij(int i, int j) const { return molecules_[i]->get_radj(j); }
  Pt get_posij(int i, int j)               {return molecules_[i]->get_posj(j);}
  Pt get_posijreal(int i, int j)
  {return molecules_[i]->get_posj_realspace(j);}
  
  shared_ptr<BaseMolecule> get_moli(int i)    { return molecules_[i]; }
  const string get_typei(int i) const {return molecules_[i]->get_move_type();}
  const double get_droti(int i) const {return molecules_[i]->get_drot();}
  const double get_dtransi(int i) const {return molecules_[i]->get_dtrans();}
  
  const int get_mol_global_idx(int type, int ty_idx)
  {
    vector<int> keys = {type, ty_idx};
    return typeIdxToIdx_[keys];
  }
   
  // Compute cutoff for force calcs
  void compute_cutoff();
  
  // Set time of simulation as what is input
  void set_time(double val) { t_ = val; }
  
  // translate every charge in molecule i by the vector dr
  void translate_mol(int i, Pt dr) { molecules_[i]->translate(dr, boxLength_); }
  
  Pt get_unwrapped_center(int i) const
      { return molecules_[i]->get_unwrapped_center(); }
  
  // rotate every charge in Molecule i
  void rotate_mol(int i, Quat qrot) { molecules_[i]->rotate(qrot); }
  void rotate_mol(int i, MyMatrix<double> rotmat)
  { return molecules_[i]->rotate(rotmat); }
  
  // Check to determine if any MoleculeSAMs are overlapping
  void check_for_overlap();
  
//  // attempt dynamics move and return whether there is a collisionv
//  virtual try_translate(
  
  // get distance vector between any two points taking into account periodic
  // boundary conditions
  Pt get_pbc_dist_vec_base(Pt p1, Pt p2);
  
  // given a distance vector, determine whether it is in the cutoff
  bool less_than_cutoff(Pt v);
  
  // write current system to PQR file, mid=-1 is print all molecules,
  // else only print one
  virtual void write_to_pqr( string outfile, int mid = -1 ) = 0;
  
  // write current system configuration to XYZ file
  virtual void write_to_xyz(ofstream &xyz_out) = 0;
};

/*
 Exception thrown when two MoleculeAMs in the system are overlapping
 */
class OverlappingMoleculeException: public exception
{
protected:
  int idx1_;
  int idx2_;
  
public:
  OverlappingMoleculeException(int idx1, int idx2)
  :idx1_(idx1), idx2_(idx2)
  {
  }
  
  virtual const char* what() const throw()
  {
    string ss;
    ss = "Molecule " + to_string(idx1_)+" & " + to_string(idx2_) + " overlap";
    return ss.c_str();
  }
};

#endif /* BaseSys_h */

//
//  System_hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef System_hpp
#define System_hpp

#include <stdio.h>
#include <vector>
#include "Constants.h"
#include "readutil.h"

using namespace std;


/*
 Class for storing relevant data about each molecule
 */
class Molecule
{
protected:
  string              type_;
  double              drot_;  // rotational diffusion coefficient
  double              dtrans_; // translational diffusion coefficients
  int                 M_;  // number of charges in this molecule
  double              a_;  // radius of this molecule
  Pt                  center_;
  vector<double>      qs_;  // magnitude of each charge in the molecule
  vector<Pt>          pos_;  // position of each charge in the molecule
  vector<double>      vdwr_; // van der waal radius of each charge
  
  
  // calculate the center of the molecule
  void calc_center();
  
  // calculate the radius of the molecule
  void calc_a();
  
  // reposition charges wrt the center
  void reposition_charges();
    
public:
  
  Molecule() { }

  // user specified radius and center
  Molecule(string type, double a, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, Pt cen, double drot_=0, double dtrans=0);
  
  // user specified radius
  Molecule(string type, double a, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, double drot_=0, double dtrans=0);
  
  // user specified center
  Molecule(string type, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, Pt cen, double drot_=0, double dtrans=0);
  
  // neither the center or radius are specified
  Molecule(string type, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, double drot_=0, double dtrans=0);
  
  const int get_m() const               { return M_; }
  const double get_a() const            { return a_; }
  const double get_qj(int j) const      { return qs_[j]; }
  Pt get_posj(int j) const              { return pos_[j]; }
  Pt get_posj_realspace(int j)          { return center_ + pos_[j]; }
  Pt get_center() const                 { return center_; }
  
  string get_type() const               { return type_; }
  
  double get_drot() const               { return drot_; }
  double get_dtrans() const             { return dtrans_; }
  
  void translate(Pt dr);
  void rotate(Quat qrot);
};


/*
 Class containing all of the relevant information for a particular system
 */
class System
{
protected:
    
  int                          N_; // number of molecules
  double                       lambda_; // average molecular radius
  vector<Molecule>             molecules_;
  Constants                    consts_;  // Constants for this system
  
  double                       boxLength_;
  double                       cutoff_;
  
  double t_;  // time in a BD simulation
  
  const double calc_average_radius() const;
  
public:
  System() { }
  
  System(Constants consts, const vector<Molecule>& mols,
         double cutoff=Constants::FORCE_CUTOFF,
         double boxlength=Constants::MAX_DIST);
  System(Constants consts, Setup setup, double cutoff=Constants::FORCE_CUTOFF);
    
  const Constants& get_consts() const      {return consts_;}
  const int get_n() const                  {return N_;}
  const double get_ai(int i) const         {return molecules_[i].get_a();}
  const double get_Mi(int i) const         {return molecules_[i].get_m();}
  const double get_qij(int i, int j) const {return molecules_[i].get_qj(j);}
  Pt get_posij(int i, int j)               {return molecules_[i].get_posj(j);}
  Molecule get_molecule(int i) const       {return molecules_[i];}
  Pt get_centeri(int i) const              {return molecules_[i].get_center();}
  double get_radi(int i) const             {return molecules_[i].get_a();}
  const double get_lambda() const          {return lambda_;}
  const string get_typei(int i) const      {return molecules_[i].get_type();}
  const double get_droti(int i) const      {return molecules_[i].get_drot();}
  const double get_dtransi(int i) const    {return molecules_[i].get_dtrans();}
  const double get_boxlength() const       {return boxLength_;}
  const double get_cutoff() const          {return cutoff_;}
  const double get_time() const            {return t_;}
  
  // Compute cutoff for force calcs
  void compute_cutoff();
  
  // Set time of simulation as what is input
  void set_time(double val) { t_ = val; }
  
  // translate every charge in molecule i by the vector dr
  void translate_mol(int i, Pt dr) { molecules_[i].translate(dr); }
  
  // rotate every charge in molecule i
  void rotate_mol(int i, Quat qrot) { molecules_[i].rotate(qrot); }
  
  // Check to determine if any molecules are overlapping
  void check_for_overlap();
  
  // get the distance vector (Point object) between two molecules, taking into
  // account periodic boundary conditions by returning the distance vector
  // between the closest image of the molecule
  Pt get_pbc_dist_vec(int i, int j);
  
  // given a distance vector, determine whether it is in the cutoff
  bool less_than_cutoff(Pt v);
  
  // return a copy of this system with a smaller set of molecules
  System get_subsystem(const vector<int> mol_idx);
  
};

/*
 Exception thrown when two molecules in the system are overlapping
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

#endif /* Setup_hpp */

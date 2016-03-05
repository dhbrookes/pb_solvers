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
#include "util.h"
#include "Constants.h"

using namespace std;


/*
 Class for storing relevant data about each molecule
 */
class Molecule
{
protected:
  
  int                 M_;  // number of charges in this molecule
  double              a_;  // radii of this molecule
  Pt                  center_;
  vector<double>      qs_;  // magnitude of each charge in the molecule
  vector<Pt>          pos_;  // position of each charge in the molecule
    
public:
  Molecule(int M, double a, vector<double> qs, vector<Pt> pos);
  Molecule(int M, vector<double> qs, vector<Pt> pos, vector<double> vwD);
  Molecule(int M, double a, vector<double> qs, vector<Pt> pos, Pt cen);
    
  const int get_m() const               { return M_; }
  const double get_a() const            { return a_; }
  const double get_qj(int j) const      { return qs_[j]; }
  Pt get_posj(int j) const              { return pos_[j]; }
  Pt get_center() const                 { return center_; }
  
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
  
  const double calc_average_radius() const;
  
public:
  System(Constants consts, const vector<Molecule>& mols);
    
  const Constants& get_consts() const { return consts_; }
  const int get_n() const { return N_; }
  const double get_ai(int i) const { return molecules_[i].get_a(); }
  const double get_Mi(int i) const { return molecules_[i].get_m(); }
  const double get_qij(int i, int j) const { return molecules_[i].get_qj(j); }
  Pt get_posij(int i, int j) { return molecules_[i].get_posj(j); }
  Molecule get_molecule(int i) const { return molecules_[i]; }
  Pt get_centeri(int i) { return molecules_[i].get_center(); }
  double get_radi(int i) { return molecules_[i].get_a(); }
  const double get_lambda()  { return lambda_; }
  
  // translate every charge in molecule i by the vector dr
  void translate_mol(int i, Pt dr) { molecules_[i].translate(dr); }
  
  // rotate every charge in molecule i
  void rotate_mol(int i, Quat qrot) { molecules_[i].rotate(qrot); }
  
  /*
   Check to determine if any molecules are overlapping
   */
  void check_for_overlap();
  
};

/*
 Exception thrown when a user-input center and radius does not encompass all
 the charges
 */
class BadCenterException: public exception
{
protected:
  double x_, y_, z_;
  double radius_;
  
public:
  BadCenterException(Pt center, double radius)
  :x_(center.x()), y_(center.y()), z_(center.z())
  {
  }
  
  virtual const char* what() const throw()
  {
    ostringstream ss;
    char buffer [50];
    printf(buffer, "Center : (%f, %f, %f) with radius : %f do \
           not encompass all charges", x_, y_, z_, radius_);
    ss << buffer << endl;
    return ss.str().c_str();
  }
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
    ostringstream ss;
    char buffer [50];
    printf(buffer, "Molecule %i and %i are overlapping", idx1_, idx2_);
    ss << buffer << endl;
    return ss.str().c_str();
  }
};

#endif /* Setup_hpp */

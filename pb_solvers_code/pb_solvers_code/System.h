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
  vector<double>      qs_;  // magnitude of each charge in the molecule
  vector<EPt>         pos_;  // position of each charge in the molecule\
    
public:
  Molecule(int M, int a, vector<double> qs, vector<EPt> pos);
    
  const int get_m() const               { return M_; }
  const double get_a() const            { return a_; }
  const double get_qj(int j) const      { return qs_[j]; }
  const EPt get_posj(int j) const       { return pos_[j]; }
  const ShPt get_sph_posj(int j) const  { return pos_[j].
                                                convert_to_spherical(); }
  
};


/*
 Class containing all of the relevant information for a particular system
 */
class System
{
protected:
    
  int                          N_; // number of molecules
  vector<Molecule>             molecules_;
  Constants                    consts_;  // Constants for this system
    
public:
  System(Constants consts, const vector<Molecule>& mols);
    
  const Constants& get_consts() const { return consts_; }
  const int get_n() const { return N_; }
  const double get_ai(int i) const { return molecules_[i].get_a(); }
  const double get_Mi(int i) const { return molecules_[i].get_m(); }
  const double get_qij(int i, int j) const { return molecules_[i].get_qj(j); }
  const EPt get_posij(int i, int j) const { return molecules_[i].get_posj(j); }
  const ShPt get_sph_posij(int i, int j) const { return molecules_[i].
                                                          get_sph_posj(j); }
  const Molecule get_molecule(int i) const { return molecules_[i]; }
  
};




#endif /* Setup_hpp */

//
//  SystemAM_h
//  pb_solvers_code
//
/*
 Copyright (c) 2015, Teresa Head-Gordon, Lisa Felberg, Enghui Yap, David Brookes
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of UC Berkeley nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SystemAM_h
#define SystemAM_h

#include <map>
#include "Constants.h"
#include "BaseSys.h"

using namespace std;


/*
 Class for storing relevant data about each MoleculeAM
 */
class MoleculeAM : public BaseMolecule
{
protected:
  Pt                  unwrappedCenter_; // unwrapped center to check for term

  vector<int>         interPol_; // List of other mols that are within 10A
  vector<int>         interAct_; // For mol, list of other are btw cutoff & 10A
  
  // calculate the center of the MoleculeAM
  Pt calc_center();
  
  // calculate the radius of the MoleculeAM
  double calc_a();
  
  // reposition charges wrt the center
  void reposition_charges();
  
public:
  
  MoleculeAM() {}

  // user specified radius and center
  MoleculeAM(string movetype, double a, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, Pt cen, int type, int typeIdx,
           double drot_=0, double dtrans=0);
  
  // user specified radius
  MoleculeAM(string movetype, double a, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, int type, int typeIdx,
           double drot_=0, double dtrans=0);
  
  // user specified center
  MoleculeAM(string movetype, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, Pt cen, int type, int typeIdx,
           double drot_=0, double dtrans=0);
  
  // neither the center or radius are specified
  MoleculeAM(string movetype, vector<double> qs, vector<Pt> pos,
           vector<double> vdwr, int type, int typeIdx,
           double drot_=0, double dtrans=0);
  
  const double get_a() const            { return as_[0]; }
  Pt get_center() const                 { return centers_[0]; }
  Pt get_unwrapped_center() const       { return unwrappedCenter_; }
  
  Pt get_cen_j(int j)                   { return centers_[0]; }

  vector<int> get_pol()                 { return interPol_;}
  vector<int> get_act()                 { return interAct_;}

  void clear_inter_pol()                { interPol_.clear(); }
  void clear_inter_act()                { interAct_.clear(); }
  void add_J_to_pol(int J)              { interPol_.push_back(J);}
  void add_J_to_interact(int J)         { interAct_.push_back(J);}
  
  bool is_J_in_pol( int J )
  {
    int j;
    for (j = 0; j < interPol_.size(); j++)
      if ( interPol_[j] == J) return true;
    
    return false;
  }
  bool is_J_in_interact( int J )
  {
    int j;
    for (j = 0; j < interAct_.size(); j++)
      if ( interAct_[j] == J) return true;
    return false;
  }  
  
  virtual void translate(Pt dr, double boxlen);
  virtual void rotate(Quat qrot);
  virtual void rotate(MyMatrix<double> rotmat);
};


/*
 Class containing all of the relevant information for a particular system
 */
class SystemAM : public BaseSystem
{
protected:
  
public:
  SystemAM() { }
  
  SystemAM(vector<shared_ptr<BaseMolecule> > mols,
         double cutoff=Constants::FORCE_CUTOFF,
         double boxlength=Constants::MAX_DIST);
   
  SystemAM(Setup setup, double cutoff=Constants::FORCE_CUTOFF);
  
  const int get_n() const                  { return N_;}
  const int get_ntype()                    { return ntype_;}
  const int get_typect(int i)              { return typect_[i];}
  const double get_ai(int i) const         { return get_aik(i, 0); }
  vector<Pt> get_allcenter() const;
  Pt get_centeri(int i) const              { return get_centerik(i, 0); }
  Pt get_unwrapped_center(int i) const
  {return molecules_[i]->get_unwrapped_center();}
  double get_radi(int i) const             { return get_aik(i, 0); }

  Pt get_pbc_dist_vec(int i, int j);

  // Interaction and polarization lists
  void add_J_to_pol_I( int i, int j) { }
  void add_J_to_interact_I(int i, int j) {}
  bool is_J_in_pol_I( int i, int j)  { return false;}
  bool is_J_in_act_I( int i, int j)  { return false;}
  vector<int> get_pol_I(int i)       { return vector<int> (1,0);}
  vector<int> get_act_I(int i)       { return vector<int> (1,0);}

  void clear_all_lists()  {  }
  
  // reset positions with input xyz file
  void reset_positions( vector<string> xyzfiles );
  
  // write current system to PQR file
  void write_to_pqr( string outfile, int mid = -1);
  
  // write current system configuration to XYZ file
  void write_to_xyz(ofstream &xyz_out);
  
};

#endif /* SystemAM_h */

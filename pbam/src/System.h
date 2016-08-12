//
//  System_hpp
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

#ifndef System_hpp
#define System_hpp

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
class System
{
protected:
    
  int                          N_; // number of molecules
  double                       lambda_; // average molecular radius
  vector<MoleculeAM>           molecules_;
  
  double                       boxLength_;
  double                       cutoff_;
  
  double t_;  // time in a BD simulation
  
  int                          ntype_;  //count of unique molecules
  vector<int>                  typect_; //count of molecules of each unique type
  map<vector<int>, int>        typeIdxToIdx_;
  
  const double calc_average_radius() const;
  
public:
  System() { }
  
  System(const vector<MoleculeAM>& mols,
         double cutoff=Constants::FORCE_CUTOFF,
         double boxlength=Constants::MAX_DIST);
   
  System(Setup setup, double cutoff=Constants::FORCE_CUTOFF);
  
  const int get_n() const                  {return N_;}
  const int get_ntype()                    {return ntype_;}
  const int get_typect(int i)              {return typect_[i];}
  const double get_ai(int i) const         {return molecules_[i].get_a();}
  const double get_Mi(int i) const         {return molecules_[i].get_m();}
  const double get_qij(int i, int j) const {return molecules_[i].get_qj(j);}
  const double get_radij(int i, int j) const {return molecules_[i].get_radj(j);}
  Pt get_posij(int i, int j)               {return molecules_[i].get_posj(j);}
  Pt get_posijreal(int i, int j)
  {return molecules_[i].get_posj_realspace(j);}
  MoleculeAM get_moli(int i) const       {return molecules_[i];}
  vector<Pt> get_allcenter() const;
  Pt get_centeri(int i) const              {return molecules_[i].get_center();}
  Pt get_unwrapped_center(int i) const
  {return molecules_[i].get_unwrapped_center();}
  double get_radi(int i) const             {return molecules_[i].get_a();}
  const double get_lambda() const          {return lambda_;}
  const string get_typei(int i) const
  {return molecules_[i].get_move_type();}
  const double get_droti(int i) const      {return molecules_[i].get_drot();}
  const double get_dtransi(int i) const    {return molecules_[i].get_dtrans();}
  const double get_boxlength() const       {return boxLength_;}
  const double get_cutoff() const          {return cutoff_;}
  const double get_time() const            {return t_;}
  
  const int get_mol_global_idx(int type, int ty_idx)
  {
    vector<int> keys = {type, ty_idx};
    return typeIdxToIdx_[keys];
  }
  
  // Compute cutoff for force calcs
  void compute_cutoff();
  
  // Set time of simulation as what is input
  void set_time(double val) { t_ = val; }
  
  // translate every charge in MoleculeAM i by the vector dr
  void translate_mol(int i, Pt dr) { molecules_[i].translate(dr, boxLength_); }
  
  // rotate every charge in MoleculeAM i
  void rotate_mol(int i, Quat qrot) { molecules_[i].rotate(qrot); }
  
  // Check to determine if any MoleculeAMs are overlapping
  void check_for_overlap();
  
  // get the distance vector (Point object) between two MoleculeAMs, taking into
  // account periodic boundary conditions by returning the distance vector
  // between the closest image of the MoleculeAM
  Pt get_pbc_dist_vec(int i, int j);
  
  // get distance vector between any two points taking into account periodic
  // boundary conditions
  Pt get_pbc_dist_vec_base(Pt p1, Pt p2);

  // Interaction and polarization lists
  void add_J_to_pol_I( int i, int j) { molecules_[i].add_J_to_pol(j);}
  void add_J_to_interact_I(int i, int j) { molecules_[i].add_J_to_interact(j);}
  
  bool is_J_in_pol_I( int i, int j)  { return molecules_[i].is_J_in_pol(j);}
  bool is_J_in_act_I( int i, int j)  { return molecules_[i].is_J_in_interact(j);}
  
  vector<int> get_pol_I(int i)       { return molecules_[i].get_pol();}
  vector<int> get_act_I(int i)       { return molecules_[i].get_act();}

  void clear_all_lists()
  {
    for ( int i = 0; i < N_; i++)
    {
      molecules_[i].clear_inter_pol();
      molecules_[i].clear_inter_act();
    }
  }  

  // given a distance vector, determine whether it is in the cutoff
  bool less_than_cutoff(Pt v);
  
  // reset positions with input xyz file
  void reset_positions( vector<string> xyzfiles );
  
  // write current system to PQR file
  void write_to_pqr( string outfile );
  
  // write current system configuration to XYZ file
  void write_to_xyz(ofstream &xyz_out);
  
};

///*
// Exception thrown when two MoleculeAMs in the system are overlapping
// */
//class OverlappingMoleculeException: public exception
//{
//protected:
//  int idx1_;
//  int idx2_;
//  
//public:
//  OverlappingMoleculeException(int idx1, int idx2)
//  :idx1_(idx1), idx2_(idx2)
//  {
//  }
//  
//  virtual const char* what() const throw()
//  {
//    string ss;
//    ss = "Molecule " + to_string(idx1_)+" & " + to_string(idx2_) + " overlap";
//    return ss.c_str();
//  }
//};

#endif /* Setup_hpp */

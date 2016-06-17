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

using namespace std;


class CGSphere
{
protected:
  Pt cen_;
  double a_;
  vector<int> ch_;  // encompassed charges
  
public:
  CGSphere() { }
  CGSphere(Pt cen, double a, vector<int> ch)
  : cen_(cen), a_(a), ch_(ch) { }
  
  Pt get_center() const { return cen_; }
  double get_a() const  { return a_; }
  int get_n() const     { return (int) ch_.size(); }
  vector<int> get_ch() const { return ch_; }
};


class Molecule
{
protected:
  string              moveType_;
  int                 type_; // int index of type of molecule, 0 based
  int                 typeIdx_; // int index of mol within given type_, 0 based
  double              drot_;  // rotational diffusion coefficient
  double              dtrans_; // translational diffusion coefficients
  
  int                 Nc_;  // number of charges in this molecule
  vector<double>      qs_;  // magnitude of each charge in the molecule
  vector<Pt>          pos_;  // position of each charge in the molecule
  vector<double>      vdwr_; // van der waal radius of each chargeGoogle
  Pt                  cog_; // Molecule center of geometry
  
  int                  Ns_;  // number of coarse grained spheres
  vector<Pt>           centers_; //coarse-grained sphere centers
  vector<double>       as_; // coarse-grained sphere radii
  vector<vector<int> > cgNeighs_; // list of indices of CG centers that neighbor
                                   // each coarse grained sphere
  vector<vector<Pt> >  cgGridPts_; // grid points on the surface of
                                  // each coarse grained sphere
  vector<vector<int> > cgGdPtExp_; // indices of grid points on the surface of
                                   // each CG sphere that are solvent exposed
  vector<vector<int> > cgGdPtBur_; // indices of grid points on the surface of
                                  // each CG sphere that are buried
  vector<vector<int> > cgCharges_; // indices of charges within each
                                   // coarse grained sphere
  vector<vector<int> > cgChargesIn_; // indices of charges within each
                                   // coarse grained sphere
  vector<vector<int> > cgChargesOut_; // indices of charges not within each
                                   // coarse grained sphere

  
  map<int, int>        chToCG_; // maps index of charge to
                                //index of its coarse-grained sphere
  
  int find_closest_center(Pt pos);
  
  /*
   Function that chooses centers for charges in the molecule
   */
  void find_centers(vector<Pt> sp, vector<Pt> np,
                    double tol_sp, int n_trials,
                    int max_trials, double beta);
  
  /*
   Perform a monte carlo search for the best center to encompass the 
   remaining chargesg
   */
  CGSphere find_best_center(vector<Pt> sp,vector<Pt> np,
                            vector<int> unbounded,
                            double tol_sp, int iter=1200, double beta=2.0);

  /* Ensure that all the CG spheres are touching */
  void check_connect();

  /* For CG sphere i, find a vector of CG spheres it's in contact w */
  vector <int> find_neighbors( int i);
  
public:
  
  Molecule() { }
  
  Molecule(const Molecule& mol);
  
  // user specified centers
  Molecule(int type, int type_idx, string movetype, vector<double> qs,
             vector<Pt> pos, vector<double> vdwr, vector<Pt> cens,
             vector<double> as, double drot=0, double dtrans=0);
  
  Molecule(int type, int type_idx, string movetype, vector<double> qs,
             vector<Pt> pos, vector<double> vdwr, vector<Pt> msms_sp,
             vector<Pt> msms_np, double tol_sp, int n_trials=10,
             int max_trials=40, double beta=2.0, double drot=0,
             double dtrans=0);
  
  void map_repos_charges();
  
  void set_type_idx(int typeidx) { typeIdx_ = typeidx; }
  void set_gridj(int j, vector<Pt> grid) {cgGridPts_[j] = grid;}
  void set_gridexpj(int j, vector<int> grid_exp) {cgGdPtExp_[j] = grid_exp;}
  void set_gridburj(int j, vector<int> grid_bur) {cgGdPtBur_[j] = grid_bur;}
  
  void translate(Pt dr, double boxlen);
  void rotate(Quat qrot); // This will rotate with respect to origin!
  void rotate(MyMatrix<double> rotmat);
  
  void calc_cog();
  
  string get_move_type() const        { return moveType_; }
  int get_type() const                { return type_; }
  int get_type_idx() const            { return typeIdx_; }
  int get_nc() const                  { return Nc_; }
  int get_ns() const                  { return Ns_; }
  int get_nc_k(int k) const           { return (int) cgCharges_[k].size(); }
  vector<int> get_neighj(int j) const { return cgNeighs_[j]; }
  vector<Pt> get_gridj(int j) const   { return cgGridPts_[j]; }
  Pt get_gridjh(int j, int h) const   { return cgGridPts_[j][h]; }
  vector<int> get_gdpt_expj(int j) const { return cgGdPtExp_[j]; }
  vector<int> get_gdpt_burj(int j) const { return cgGdPtBur_[j]; }
  
  vector<int> get_ch_allin_k(int k)   { return cgChargesIn_[k]; }
  vector<int> get_ch_allout_k(int k)  { return cgChargesOut_[k]; }
  
  int get_ch_k_alpha(int k, int alpha){ return cgCharges_[k][alpha]; }
  double get_drot() const             { return drot_; }
  double get_dtrans() const           { return dtrans_; }
  Pt get_posj(int j) const            { return pos_[j]; }
  
  Pt get_posj_realspace(int j)        { return pos_[j] + centers_[chToCG_[j]];}
  Pt get_centerk(int k) const         { return centers_[k]; }
  Pt get_cog() const                  { return cog_;}
  const double get_qj(int j) const    { return qs_[j]; }
  const double get_radj(int j) const  { return vdwr_[j]; }
  const double get_ak(int k) const    { return as_[k]; }
  const int get_cg_of_ch(int j)       { return chToCG_[j]; }
  
  /* Choose a random orientation for a Pt vector  */
  Pt random_pt();
  
  // random number from normal distribution
  double random_norm();
  
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
  
  double                       boxLength_;
  double                       cutoff_;
  
  double t_;  // time in a BD simulation
  
  int                          ntype_;  //count of unique molecule types
  vector<int>                  typect_; //count of molecule of each unique type
  map<vector<int>, int>        typeIdxToIdx_;
  
  const double calc_average_radius() const;
  
//  Molecule build_type_mol(int type, Setup setup);
  
public:
  System() { }
  
  System(const vector<Molecule>& mols,
         double cutoff=Constants::FORCE_CUTOFF,
         double boxlength=Constants::MAX_DIST);
  
  System(Setup setup, double cutoff=Constants::FORCE_CUTOFF);
  
  const int get_n() const                  {return N_;}
  const int get_ntype()                    {return ntype_;}
  const int get_typect(int i)              {return typect_[i];}
  const double get_aik(int i, int k) const {return molecules_[i].get_ak(k);}
  const double get_Nc_i(int i) const       {return molecules_[i].get_nc();}
  const double get_Ns_i(int i) const       {return molecules_[i].get_ns();}
  const double get_qij(int i, int j) const {return molecules_[i].get_qj(j);}
  const double get_droti(int i) const      {return molecules_[i].get_drot();}
  const double get_dtransi(int i) const    {return molecules_[i].get_dtrans();}
  const double get_boxlength() const       {return boxLength_;}
  const double get_cutoff() const          {return cutoff_;}
  const double get_time() const            {return t_;}
  const double get_lambda() const          {return lambda_;}
  Molecule get_molecule(int i) const       {return molecules_[i];}
  Pt get_cogi(int i) const                {return molecules_[i].get_cog();}
  Pt get_posij(int i, int j)               {return molecules_[i].get_posj(j);}
  Pt get_centerik(int i, int k) const   {return molecules_[i].get_centerk(k);}
  const string get_typei(int i) const   {return molecules_[i].get_move_type();}
  const double get_radij(int i, int j) 
                                     const {return molecules_[i].get_radj(j);}
  Pt get_posijreal(int i, int j) {return molecules_[i].get_posj_realspace(j);}
  
  Pt get_gridijh(int i, int j, int h) const
        { return molecules_[i].get_gridjh(j, h); }
  vector<int> get_gdpt_expij(int i, int j) const
        { return molecules_[i].get_gdpt_expj(j); }
  
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
  void translate_mol(int i, Pt dr) { molecules_[i].translate(dr, boxLength_); }
  
  // rotate every charge in molecule i
  void rotate_mol(int i, Quat qrot) { molecules_[i].rotate(qrot); }
  
  // Check to determine if any molecules are overlapping
  void check_for_overlap();
  
  // get distance vector between any two points taking into account periodic
  // boundary conditions
  Pt get_pbc_dist_vec_base(Pt p1, Pt p2);
  
  // given a distance vector, determine whether it is in the cutoff
  bool less_than_cutoff(Pt v);
  
  // write current system to PQR file
  void write_to_pqr( string outfile );
  
  // write current system configuration to XYZ file
  void write_to_xyz(ofstream &xyz_out);
  
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

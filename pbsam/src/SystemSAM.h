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
#include <memory>
#include "Constants.h"
#include "BaseSys.h"

using namespace std;

namespace pbsolvers
{

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


class MoleculeSAM : public BaseMolecule
{
protected:
  Pt                  cog_; // MoleculeSAM center of geometry
  Pt                  cog_unwrapped_; // MoleculeSAM center of geometry unwrapp
  
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

  vector<vector<vector<int> > > interPol_; // For each sph in mol, list of
                                       // mol/sph pairs that are within 10A
  vector<vector<vector<int> > > interAct_; // For each sph in mol, list of
                                      // mol/sph pairs that are btw 100 & 10A
  
  map<int, int>        chToCG_; // maps index of charge to
                                //index of its coarse-grained sphere
  
  int find_closest_center(Pt pos);
  
  /*
   Function that chooses centers for charges in the MoleculeSAM
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
                            double tol_sp, double beta=2.0);

  /* Ensure that all the CG spheres are touching */
  void check_connect();

  /* For CG sphere i, find a vector of CG spheres it's in contact w */
  vector <int> find_neighbors( int i);
  
  // random number from normal distribution
  double random_norm();
  
  /* Choose a random orientation for a Pt vector  */
  Pt random_pt();

  void map_repos_charges();
  
public:
  
  MoleculeSAM() { }
  
  MoleculeSAM(const MoleculeSAM& mol);
  
  // user specified centers
  MoleculeSAM(int type, int type_idx, string movetype, vector<double> qs,
             vector<Pt> pos, vector<double> vdwr, vector<Pt> cens,
             vector<double> as, double drot=0, double dtrans=0);
  
  MoleculeSAM(int type, int type_idx, string movetype, vector<double> qs,
             vector<Pt> pos, vector<double> vdwr, vector<Pt> msms_sp,
             vector<Pt> msms_np, double tol_sp, int n_trials=10,
             int max_trials=40, double beta=2.0, double drot=0,
             double dtrans=0);
  
  MoleculeSAM(int type, int type_idx, string movetype, vector<double> qs,
             vector<Pt> pos, vector<double> vdwr, string msms_fil,
             double tol_sp, double drot=0, double dtrans = 0, int n_trials=10,
             int max_trials=40, double beta=2.0);
  
  void set_type_idx(int typeidx) { typeIdx_ = typeidx; }

  
  void add_Jl_to_interk( int k, int J, int l) {interPol_[k].push_back({J,l});}
  void add_Jl_to_inter_act_k( int k, int J, int l)
  {interAct_[k].push_back({J,l});}
  
  bool is_J_in_interk( int k, int J )
  {
    int j;
    for (j = 0; j < interPol_[k].size(); j++)
      if ( interPol_[k][j][0] == J) return true;
    
    return false;
  }
  
  void set_gridj(int j, vector<Pt> grid) {cgGridPts_[j] = grid;}
  void set_gridexpj(int j, vector<int> grid_exp) {cgGdPtExp_[j] = grid_exp;}
  void set_gridburj(int j, vector<int> grid_bur) {cgGdPtBur_[j] = grid_bur;}
  
  void translate(Pt dr, double boxlen);
  void rotate(Quat qrot); // This will rotate with respect to origin!
  void rotate(MyMatrix<double> rotmat);
  
  void calc_cog();
  void write_pqr(string outfile);

  int get_nc_k(int k) const           { return (int) cgCharges_[k].size(); }
  vector<int> get_neighj(int j) const { return cgNeighs_[j]; }
  vector<Pt> get_gridj(int j) const   { return cgGridPts_[j]; }
  Pt get_gridjh(int j, int h) const   { return cgGridPts_[j][h]; }
  vector<int> get_gdpt_expj(int j) const { return cgGdPtExp_[j]; }
  vector<int> get_gdpt_burj(int j) const { return cgGdPtBur_[j]; }
  vector<int> get_ch_allin_k(int k)   { return cgChargesIn_[k]; }
  vector<int> get_ch_allout_k(int k)  { return cgChargesOut_[k]; }
  int get_ch_k_alpha(int k, int alpha){ return cgCharges_[k][alpha]; }
  Pt get_cen_j(int j)                 { return centers_[chToCG_[j]]; }
  Pt get_cog() const                  { return cog_;}
  Pt get_unwrapped_center() const     { return cog_unwrapped_; }
  const int get_cg_of_ch(int j)       { return chToCG_[j]; }
  vector<vector<int> > get_inter_act_k(int k) {return interAct_[k]; }

};

/*
 Class containing all of the relevant information for a particular system
 */
class SystemSAM : public BaseSystem
{
protected:
  vector<vector<double> >      min_dist_; // minimum dist between mols
  
public:
  SystemSAM() { }
  
  SystemSAM(vector<shared_ptr<BaseMolecule> > mols,
         double cutoff=Constants::FORCE_CUTOFF,
         double boxlength=Constants::MAX_DIST);
  
  SystemSAM(Setup setup, double cutoff=Constants::FORCE_CUTOFF);
  
  vector<int> get_all_Ik()
  {
    vector<int> ns_i;
    for ( int i = 0; i< N_; i++) ns_i.push_back(get_Ns_i(i));
    return ns_i;
  }
  
  Pt get_cogi(int i) const                {return molecules_[i]->get_cog();}

  const double get_radij(int i, int j)
                                     const {return molecules_[i]->get_radj(j);}
  Pt get_posijreal(int i, int j) {return molecules_[i]->get_posj_realspace(j);}
  
  Pt get_unwrapped_center(int i) const
  {return molecules_[i]->get_unwrapped_center();}
  
  Pt get_gridijh(int i, int j, int h) const
        { return molecules_[i]->get_gridjh(j, h); }
  vector<int> get_gdpt_expij(int i, int j) const
        { return molecules_[i]->get_gdpt_expj(j); }
  
  const int get_mol_global_idx(int type, int ty_idx)
  {
    vector<int> keys = {type, ty_idx};
    return typeIdxToIdx_[keys];
  }
  
  double get_min_dist(int I, int J) { return min_dist_[min(I, J)][max(I, J)];}
  
  // Compute min sphere center 2 center distance between 2 mols
  double calc_min_dist(int I, int J);
  
  // Save min distances for all MoleculeSAM pairs
  void save_min_dist();
  
  //Copy surface integral points from molecule i to molecule j. MUST
  // be the same type!
  void copy_grid( int i, int j)
  {
    for (int k = 0; k < molecules_[i]->get_ns(); k++)
    {
      molecules_[j]->set_gridj(k, molecules_[i]->get_gridj(k));
      molecules_[j]->set_gridburj(k, molecules_[i]->get_gdpt_burj(k));
      molecules_[j]->set_gridexpj(k, molecules_[i]->get_gdpt_expj(k));
    }
  }
  
  // Reset positions for new BD trajectory
  void reset_positions(vector<string> xyzfiles);
  
  // write current system to PQR file, mid=-1 is print all MoleculeSAMs,
  // else only print one
  void write_to_pqr(string outfile, int mid = -1 );
  
  // write current system configuration to XYZ file
  void write_to_xyz(ofstream &xyz_out);
  
};

/*
 Exception thrown when two molecules in the system are overlapping
 */
class OverlappingMoleculeSAMException: public exception
{
protected:
  int idx1_;
  int idx2_;
  mutable string ss_;
  
public:
  OverlappingMoleculeSAMException(int idx1, int idx2)
  :idx1_(idx1), idx2_(idx2), ss_()
  {
  }
  
  virtual const char* what() const throw()
  {
    ss_ = "MoleculeSAM " + to_string(idx1_)+" & " + to_string(idx2_) + " overlap";
    return ss_.c_str();
  }
};

} /* namespace pbsolvers */
#endif /* Setup_hpp */

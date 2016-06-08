//
// setup.h
// pbam
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

#ifndef SETUP_H
#define SETUP_H

#include "readutil.h"

using namespace std;

class Setup
{
protected:
  
  static const double MAX_DIST;  // maximum distance for cutoff, box length
  
  string  units_; // the units desired for output
  int     ompThreads_;
  int     nType_;  		// Number of different molecule types
  int     PBCs_;			// PBC in 0, pseudo-2, or 3 dimensions
  int     maxtime_;
  int     ntraj_;
  int     srand_;			// random seed
  double  saltConc_;
  double  blen_; 		// boxlength for PBC
  double  idiel_;
  double  sdiel_;  // dielectric constant win molecule and of solvent
  double  iKbT_;
  double  temp_;
  double  kappa_;
  bool    orientRand_; // flag for creating random orientations for mols

  // make spheres settings:
  double    sphBeta_;
  double    tolSP_;
  int       nTrials_;
  int       maxTrials_;
  
  // for electrostatics runtype
  int             gridPts_;  // number of voxels to compute for each dim
  int             gridCt_;  // number of grid files to write
  vector<string>  axis_;    // For grid print, axis desired
  vector<double>  axLoc_;   // Location along given axis
  vector< string> potOutfnames_; // Vector of outfiles, [0] = dx, rest = grid
  
  // for dynamics runs
  int                   numTerm_;   //number of termination conditions
  vector<string>        termtype_;  // type of each term ('time', 'x', 'y', 'z',
                                    //'r' or 'contact')
  vector<string>        confiles_;  // contact files for contact termination
                                    //conditions
  vector<double>        conpads_;   // pads for contact termination conditions
  vector<double>        termvals_;  // value for each termination condition
  vector<vector<int> >  termmols_;  // vector of molecule ids
  bool                  andCombine_;  //if true, term conds will
                                      //combine w/ 'and', otherwise 'or'
  
  vector<int>             nTypenCount_;  // number of mols of each type
  vector<vector<double> > typeDiff_;  // Dtr, Drot each type, size [Ntype][2]
  vector<string>          typeDef_;  // For each type, type is stat, rot or move
  vector<string>          runSpecs_;  //include run type [0] (electrost/bd)
                                      //and runname [1]
  vector<string>          pqr_names_;  // PQR file names
  vector<string>          surfNames_;  // MSMS surface file names
  vector<vector<string> > xyz_names_;  // XYZ file names
  vector<vector<bool> >   isTransRot_;
  
  vector<string> mbdfile_loc_; // location of names for manybd data output
  

  // input file reading methods:
  void read_infile(string infile);
  vector<string> split(string str, char delim);
  void findLines(string fline);
  void findKeyword(vector<string> fline);
  void resizeVecs();
  
  // basic settings:
  void setRunType( string runt )      { runSpecs_[0] = runt; }
  void setRunName( string runn )      { runSpecs_[1] = runn; }
  void setUnits( string units )       { units_ = units;}
  void setIDiel( double idiel )       { idiel_ = idiel; }
  void setSDiel( double sdiel )       { sdiel_ = sdiel;}
  void setTemp( double temp )         { temp_ = temp;}
  void setRand( int rand )            { srand_ = rand; }
  void setRandOrient()                { orientRand_ = true; }
  void setOMP( int ompT )             { ompThreads_ = ompT ; }
  void setSaltCon( double saltCon )   { saltConc_ = saltCon; }
  void setNType( int numType )        { nType_ = numType; }
  void setPBCT( int pbc )             { PBCs_ = pbc; }
  void setBoxl( double boxl )         { blen_ = boxl; }
  void setMaxTime( int maxt )         { maxtime_ = maxt; }
  void setKappa( double kappa )       { kappa_ = kappa; }
  void set_tol_sp(double tolsp)       { tolSP_ = tolsp; }
  void set_sph_beta(double sphbeta)   { sphBeta_ = sphbeta; }
  void set_n_trials(int n)            { nTrials_ = n; }
  void set_max_trials(int n)          { maxTrials_ = n; }
  
  // three body settings:
  void set2BDLoc( string fileloc )    { mbdfile_loc_[0] = fileloc;}
  void set3BDLoc( string fileloc )    { mbdfile_loc_[1] = fileloc;}
  
  // electrostatics settings:
  void setDXoutName( string dx)             { potOutfnames_[0] = dx;}
  void set_3dmap_name( string ht)           { potOutfnames_[1] = ht;}
  void setGridOutName( int i, string grid)  { potOutfnames_[i+1] = grid;}
  void setGridPts( int gridP )              { gridPts_ = gridP; }
  void setGridCt( int gridC )               { gridCt_ = gridC; }
  void setGridAx( int i, string ax)         { axis_[i-1] = ax;}
  void setGridAxLoc( int i, double axLoc)   { axLoc_[i-1] = axLoc;}
  
  
  //dynamics settings:
  void set_numterms(int n) { numTerm_ = n; }
  void setNTraj( int ntraj )          { ntraj_ = ntraj; }
  void resize_termcond(int n)
  {
    termtype_.resize(n);
    termmols_.resize(n);
    for ( int i = 0; i < n; i++) termmols_[i].resize(2);
    termvals_.resize(n);
  }
  
  void add_termcond(int i, string type, vector<int> mol_idx, double val)
  {
    termtype_[i] = type;
    for (int n = 0; n < 2; n++) termmols_[i][n] = mol_idx[n];
    termvals_[i] = val;
  }
  
  void set_term_combine(string type)
  {
    if (type=="and") andCombine_ = true;
    else if (type=="or") andCombine_ = false;
    else andCombine_ = false;
  }
  
  // file paths:
  void setTypeNCount( int typeCount, int count )
  { nTypenCount_[typeCount] = count; }
  
  void setTypeNDef( int typeCount, string definit )
  { typeDef_[typeCount] = definit; }
  
  void setTypeNDtr( int typeCount, double dTR )
  { typeDiff_[typeCount][0] = dTR; }
  
  void setTypeNDrot( int typeCount, double dRot )
  { typeDiff_[typeCount][1] = dRot; }
  
  void setTypeNPQR( int typeCount, string pqr )
  { pqr_names_[typeCount] = pqr; }
  
  void setTypeNSurf( int typeCount, string path )
  { surfNames_[typeCount] = path; }
  
  void setTypeNXYZ( int typeCount, int traj, string xyz )
  { xyz_names_[typeCount][traj] = xyz; }
  
  void setTypeNisTransRot( int typeCount, int traj, bool istr )
  { isTransRot_[typeCount][traj] = istr; }
  
public:
  Setup(string infile);

  // Not pretty, but getting necessary inputs from APBS
  Setup(double temp, double salt_conc, double int_diel, double solv_diel, 
        int nmol, string runtype, string runname, bool randorient, double boxl,
        int pbc_type, int gridpts, string map3d, int g2dct, 
        vector<string> grid2Dfn, vector <string> grid2Dax, 
        vector<double> grid2Dloc, string dxnam, int ntraj, bool termcomb, 
        vector<string> difftype, vector<vector<double> > diffcon,
        vector<string> termcond, vector<double> termval, 
        vector<vector <int > > termnu, vector<string> confil,
        vector<double> conpad, vector<vector <string> > xyzfil);
  
  // electrostatics
  string getDXoutName(  )         { return potOutfnames_[0];}
  string get_3dmap_name( )        { return potOutfnames_[1];}
  string getGridOutName( int i )  { return potOutfnames_[i+2];}
  
  int getGridPts() { return gridPts_; }
  int getGridCt() { return gridCt_; }
  string getGridAx( int i ) { return axis_[i];}
  double getGridAxLoc( int i ) { return axLoc_[i];}
  
  // Dynamics termination conds
  int get_numterms( )                 { return numTerm_; }
  string get_termtype( int i)         { return termtype_[i]; }
  vector<int> get_termMolIDX( int i)  { return termmols_[i]; }
  double get_termval( int i)          { return termvals_[i];}
  string get_confile(int j)           { return confiles_[j]; }
  double get_conpad(int j)            { return conpads_[j]; }
  bool get_andCombine( )              { return andCombine_; }
  int getMaxTime()                    { return maxtime_; }
  int getNTraj()                      { return ntraj_; }
  
  // threebody
  string get2BDLoc()               { return mbdfile_loc_[0]; }
  string get3BDLoc()               { return mbdfile_loc_[1]; }
  vector<string> getMBDLoc()       { return mbdfile_loc_; }
  
  // retreive basic settings:
  string getRunType()              { return runSpecs_[0]; }
  string getRunName()              { return runSpecs_[1]; }
  string getUnits()                { return units_; }
  string getTypeNXYZ(int type)     { return xyz_names_[type][0]; }
  bool get_randOrient()            { return orientRand_; }
  int getThreads()                 { return ompThreads_; }
  int getNType()                   { return nType_; }
  int getTypeNCount(int type)      { return nTypenCount_[type]; }
  int getPBCs()                    { return PBCs_; }
  int get_n_trials()               { return nTrials_; }
  int get_max_trials()             { return maxTrials_; }
  double getBLen()                 { return blen_; }
  double getIDiel()                { return idiel_; }
  double getSDiel()                { return sdiel_; }
  double getSaltConc()             { return saltConc_; }
  double getTemp()                 { return temp_; }
  double getDtr( int n )           { return typeDiff_[n][0]; }
  double getDrot( int n )          { return typeDiff_[n][1]; }
  double getKappa()                { return kappa_; }
  double getIKbT()                 { return iKbT_; }
  double get_tol_sp()              { return tolSP_; }
  double get_sph_beta ()           { return sphBeta_; }
  vector<int> get_type_nct()       { return nTypenCount_;}

  // retrieve files:
  string getTypeNDef(int type)                { return typeDef_[type]; }
  string getTypeNPQR(int type)                { return pqr_names_[type]; }
  string getTypeNSurf(int type)               { return surfNames_[type]; }
  string getTypeNXYZ(int type, int traj)      { return xyz_names_[type][traj];}
  bool getTypeIsTransRot(int type, int traj)  { return isTransRot_[type][traj];}
  bool getTypeIsTransRot(int type)            { return isTransRot_[type][0]; }
  vector<string> get_trajn_xyz(int traj)
  {
    vector<string> traj_xyz;
    for (int i = 0; i < nType_; i ++) traj_xyz.push_back(xyz_names_[i][traj]);
    return traj_xyz;
  }
  
  // check input arguments and throw BadInputException if they dont check out
  void check_inputs();
  
};

class BadInputException: public exception
{
protected:
  vector<string> problems_;
  
public:
  BadInputException(vector<string> problems)
  :problems_(problems)
  {
  }
  
  virtual const char* what() const throw()
  {
    string ss;
    ss = "The following problems were found in your input file:\n";
    for (int i = 0; i < problems_.size(); i++)
    {
      ss += to_string(i) + ": " + problems_[i] + "\n";
    }
    return ss.c_str();
  }
};

#endif

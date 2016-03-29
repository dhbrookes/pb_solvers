
#ifndef SETUP_H
#define SETUP_H

#include "readutil.h"

using namespace std;

class Setup
{
protected:
  
  static const double MAX_DIST;  // maximum distance for cutoff, box length
  
  int ompThreads_;
  double saltConc_;
  int nType_;  		// Number of different molecule types
  int PBCs_;			// PBC in 0, pseudo-2, or 3 dimensions
  double blen_; 		// boxlength for PBC
  int maxtime_;
  int ntraj_;
  
  // for electrostatics runtype
  vector< string> potOutfnames_; // Vector of outfiles, [0] = dx, rest = grid
  
  int gridPts_; // number of voxels to compute for each dim
  int gridCt_; // number of grid files to write
  vector<string> axis_;  // For grid print, axis desired
  vector<double> axLoc_; // Location along given axis
  
  int numTerm_;  //number of termination conditions
  vector<string> termtype_; // type of each term ('time', 'x', 'y', 'z' or 'r')
  vector<vector<int> > termmols_; // vector of molecule ids
  vector<double> termvals_; // value for each termination condition
  bool andCombine_;  //if true, termination conditions will combine with 'and', otherwise with 'or'
  
  double idiel_;
  double sdiel_;  // dielectric constant win molecule and of solvent
  double temp_;
  int srand_;			// random seed
  bool orientRand_; // flag for creating random orientations for mols
  
  vector<int> nTypenCount_; // Array for each of mol types, how many mols
  vector<vector<double> > typeDiff_; // Dtr, Drot each type, size [Ntype][2]
  vector<string> typeDef_; 		// For each type, type is stat, rot or move
  vector<string> runSpecs_;	//include run type [0] (electrost/bd) & runname [1]
  vector<vector<string> > molfnames_;  // file names
  
  vector<string> mbdfile_loc_; // location of names for manybd data output
  
  double kappa_;
  double iKbT_;
  
  string units_; // the units desired for output
  
  // input file reading methods:
  void read_infile(string infile);
  vector<string> split(string str, char delim);
  void findLines(string fline);
  void findKeyword(vector<string> fline);
  
  //'electrostatics' or 'dynamics' (for RunType)
  void setRunType( string runt ) {runSpecs_[0] = runt;}
  void setRunName( string runn ) {runSpecs_[1] = runn;}
  void setUnits( string units )  {units_ = units;}
  void resizeVecs();
  
  // setting values for electrostatics run
  void setDXoutName( string dx) { potOutfnames_[0] = dx;}
  void setGridOutName( int i, string grid) { potOutfnames_[i] = grid;}
  
  void setGridPts( int gridP ) { gridPts_ = gridP; }
  void setGridCt( int gridC ) { gridCt_ = gridC; }
  void setGridAx( int i, string ax) { axis_[i-1] = ax;}
  void setGridAxLoc( int i, double axLoc) { axLoc_[i-1] = axLoc;}
  
  //dynamics settings
  void set_numterms(int n) { numTerm_ = n; }
  void add_termcond(string type, vector<int> mol_idx, double val)
  {
    termtype_.push_back(type);
    termmols_.push_back(mol_idx);
    termvals_.push_back(val);
  }
  
  void set_term_combine(string type)
  {
    if (type=="and") andCombine_ = true;
    else if (type=="or") andCombine_ = false;
    else andCombine_ = false;
  }
  
  // setting details for 3bd
  void set2BDLoc( string fileloc ) { mbdfile_loc_[0] = fileloc;}
  void set3BDLoc( string fileloc ) { mbdfile_loc_[1] = fileloc;}
  
  //
  void setOMP( int ompT ) { ompThreads_ = ompT ; }
  void setSaltCon( double saltCon )
  { saltConc_ = saltCon; }
  void setNType( int numType ) { nType_ = numType; }
  void setPBCT( int pbc ){ PBCs_ = pbc; }
  void setBoxl( double boxl ){ blen_ = boxl; }
  void setMaxTime( int maxt ){ maxtime_ = maxt; }
  
  void setIDiel( double idiel ) { idiel_ = idiel; }
  void setSDiel( double sdiel ) { sdiel_ = sdiel;}
  void setTemp( double temp ) { temp_ = temp;}
  void setRand( int rand ) { srand_ = rand; }
  void setRandOrient()          { orientRand_ = true; }
  void setNTraj( int ntraj ){ ntraj_ = ntraj; }
  void setKappa( double kappa ) { kappa_ = kappa; }
  
  void setTypeNCount( int typeCount, int count )
  { nTypenCount_[typeCount] = count; }
  void setTypeNDef( int typeCount, string definit )
  { typeDef_[typeCount] = definit; }
  void setTypeNDtr( int typeCount, double dTR )
  { typeDiff_[typeCount][0] = dTR; }
  void setTypeNDrot( int typeCount, double dRot )
  { typeDiff_[typeCount][1] = dRot; }
  
  void setTypeNPQR( int typeCount, string pqr )
  { molfnames_[typeCount][0] = pqr; }
  
  void setTypeNXYZ( int typeCount, string xyz )
  { molfnames_[typeCount][1] = xyz; }
  
public:
  Setup(string infile);
  
  string getRunType()              { return runSpecs_[0]; }
  string getRunName()              { return runSpecs_[1]; }
  string getUnits()                { return units_; }
  
  // electrostatics
  string getDXoutName(  ) { return potOutfnames_[0];}
  string getGridOutName( int i ) { return potOutfnames_[i+1];}
  
  int getGridPts() { return gridPts_; }
  int getGridCt() { return gridCt_; }
  string getGridAx( int i ) { return axis_[i];}
  double getGridAxLoc( int i ) { return axLoc_[i];}
  
  // Dynamics termination conds
  int get_numterms( )              { return numTerm_; }
  string get_termtype( int i)      { return termtype_[i]; }
  vector<int> get_termMolIDX( int i) { return termmols_[i]; }
  double get_termval( int i)       { return termvals_[i];}
  bool get_andCombine( )           { return andCombine_; }
  
  
  // threebody
  string get2BDLoc()               { return mbdfile_loc_[0]; }
  string get3BDLoc()               { return mbdfile_loc_[1]; }
  vector<string> getMBDLoc()       { return mbdfile_loc_; }
  
  int getThreads()                 { return ompThreads_; }
  int getNType()                   { return nType_; }
  int getPBCs()                    { return PBCs_; }
  double getBLen()                 { return blen_; }
  double getIDiel()                { return idiel_; }
  double getSDiel()                { return sdiel_; }
  double getSaltConc()             { return saltConc_; }
  double getTemp()                 { return temp_; }
  int getMaxTime()                 { return maxtime_; }
  int getNTraj()                   { return ntraj_; }
  double getDtr( int n )           { return typeDiff_[n][0]; }
  double getDrot( int n )          { return typeDiff_[n][1]; }
  int getTypeNCount(int type)      { return nTypenCount_[type]; }
  string getTypeNDef(int type)     { return typeDef_[type]; }
  string getTypeNPQR(int type)     { return molfnames_[type][0]; }
  string getTypeNXYZ(int type)     { return molfnames_[type][1]; }
  double getKappa()                { return kappa_; }
  double getIKbT()                 { return iKbT_; }
  int get_ntype()                  { return nType_; }
  
  bool get_randOrient()            { return orientRand_; }
  
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

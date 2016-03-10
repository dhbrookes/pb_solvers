
#ifndef SETUP_H
#define SETUP_H

#include "readutil.h"
#include <fstream>

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
  
  double idiel_;
  double sdiel_;  // dielectric constant win molecule and of solvent
  double temp_;
  int srand_;			// random seed
  
  vector<int> nTypenCount_; // Array for each of mol types, how many mols
  vector<vector<double> > typeDiff_; // Dtr, Drot each type, size [Ntype][2]
  vector<string> typeDef_; 		// For each type, type is stat, rot or move
  vector<string> runSpecs_;	//include run type [0] (electrost/bd) & runname [1]
  vector<vector<string> > molfnames_;  // file names
  
  double kappa_;
  double iKbT_;
  
  // input file reading methods:
  void read_infile(string infile);
  vector<string> split(string str, char delim);
  void findLines(string fline);
  void findKeyword(vector<string> fline);
  
  //'electrostat' or 'dynamics' (for RunType)
  void setRunType( string runt ) {runSpecs_[0] = runt;}
  void setRunName( string runn ) {runSpecs_[1] = runn;}
  void resizeVecs();
  
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
  
  int getThreads()                 { return ompThreads_; }
  int getType()                    { return nType_; }
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
  
};

#endif

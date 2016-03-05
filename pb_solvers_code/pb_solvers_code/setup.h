#ifndef SETUP_H
#define SETUP_H

#include <cfloat>
#include "float.h"
#include "util.h"
#include "Constants.h"
#include "System.h"
#include "readutil.h"
//#include "constant.h"

using namespace std;

class CSetup
{
public:
  CSetup();
  ~CSetup() {};

  void printSetupClass( );

  void setRunType( string runt ) {m_runSpecs[0] = runt ;}  //'electrostat' or 'dynamics'
  void setRunName( string runn ) {m_runSpecs[1] = runn ;}
  void resizeVecs();

  void setOMP( int ompT ) { m_ompThreads = ompT ; }
  void setSaltCon( double saltCon ) 
									{ m_saltConc = saltCon; }
  void setNType( int numType ) { m_nType = numType; }
  void setPBCT( int pbc ){ m_PBCs = pbc; }
  void setBoxl( double boxl ){ m_blen = boxl; }
  void setMaxTime( int maxt ){ m_maxtime = maxt; }

  void setIDiel( double idiel ) { m_idiel = idiel; } 
  void setSDiel( double sdiel ) { m_sdiel = sdiel;}
  void setTemp( double temp ) { m_temp = temp;}
  void setRand( int rand ) { m_srand = rand; }
  void setNTraj( int ntraj ){ m_ntraj = ntraj; }
  void setKappa( double kappa ) { m_kappa = kappa; }

  void setTypeNCount( int typeCount, int count ) 
              { m_nTypenCount[typeCount] = count; }
  void setTypeNDef( int typeCount, string definit ) 
              { m_typeDef[typeCount] = definit; }
  void setTypeNDtr( int typeCount, double dTR ) 
              { m_typeDiff[typeCount][0] = dTR; }
  void setTypeNDrot( int typeCount, double dRot ) 
              { m_typeDiff[typeCount][1] = dRot; }

  void setTypeNPQR( int typeCount, string pqr ) 
              { m_molfnames[typeCount][0] = pqr; }

  void setTypeNXYZ( int typeCount, string xyz ) 
              { m_molfnames[typeCount][1] = xyz; }
	
	string getRunType()              { return m_runSpecs[0]; }
	string getRunName()              { return m_runSpecs[1]; }

  int getThreads()                 { return m_ompThreads; }
	int getType()                    { return m_nType; }
	int getPBCs()                    { return m_PBCs; }
  double getBLen()                 { return m_blen; }
  double getIDiel()                { return m_idiel; }
  double getSDiel()                { return m_sdiel; }
  double getSaltConc()             { return m_saltConc; }
  double getTemp()                 { return m_temp; }
  int getMaxTime()                 { return m_maxtime; }


  int getNTraj()                   { return m_ntraj; }


  double getDtr( int n )           { return m_typeDiff[n][0]; }
  double getDrot( int n )          { return m_typeDiff[n][1]; }

	int getTypeNCount(int type)      { return m_nTypenCount[type]; }
	string getTypeNDef(int type)     { return m_typeDef[type]; }
	string getTypeNPQR(int type)     { return m_molfnames[type][0]; }

	string getTypeNXYZ(int type)     { return m_molfnames[type][1]; }

  double getKappa()                { return m_kappa; }
  double getIKbT()                 { return m_IKbT; }
  
  // setup and return a Constants object
  Constants setup_constants();
  System setup_system(Constants consts);

private:
  int m_ompThreads;
  double m_saltConc;
  int m_nType;  		// Number of different molecule types
  int m_PBCs;			// PBC in 0, pseudo-2, or 3 dimensions
  double m_blen; 		// boxlength for PBC
  int m_maxtime;

  int m_ntraj;

  double m_idiel, m_sdiel;// dielectric constant win molecule and of solvent
  double m_temp; 
  int m_srand;			// random seed

  vector<int> m_nTypenCount; // Array for each of mol types, how many mols
  vector<vector<double> > m_typeDiff; // Dtr,Drot  each type, size [Ntype][2]
  vector<string> m_typeDef; 		// For each type, type is stat, rot or move
  string m_runSpecs[2];	//include run type [0],
                        // pbsam/bd and the runname [1]
  vector<vector<string> > m_molfnames;  // file names 

  double m_kappa, m_IKbT;

}; // end CSetup

#endif

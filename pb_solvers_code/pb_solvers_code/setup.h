#ifndef SETUP_H
#define SETUP_H

#include <cfloat>
#include "float.h"
#include "util.h"
#include "Constants.h"
//#include "constant.h"

using namespace std;

class CSetup
{
public:
  CSetup(Constants consts);
  ~CSetup() {};

  void resetKappaFact2( );
  void printSetupClass( );

  void setRunType( string runt ) {m_runSpecs[0] = runt ;}  //'electrostat' or 'dynamics'
//  void setSimType( string simt ) {m_runSpecs[1] = simt ;}
  void setRunName( string runn ) {m_runSpecs[1] = runn ;}
  void resizeVecs();

  void setOMP( int ompT ) { m_ompThreads = ompT ; }
  void setSaltCon( double saltCon ) 
									{ m_saltConc = saltCon; resetKappaFact2(); }
  void setNType( int numType ) { m_nType = numType; }
  void setPBCT( int pbc ){ m_PBCs = pbc; }
  void setBoxl( double boxl ){ m_blen = boxl; }
//  void setCutoffInter( double intermolCutoff ) 
//											{ m_cutInter = intermolCutoff ; }
//  void setCutoffInteract( double interactmolCutoff ) 
//											{ m_cutInteract = interactmolCutoff ; }
//  void setCutoffIntra( double intramolCutoff ) 
//											{ m_cutIntra = intramolCutoff ; }
  void setMaxTime( int maxt ){ m_maxtime = maxt; }

//  void setDistance( double distance ){ m_distance = distance; }

  void setIDiel( double idiel ) { m_idiel = idiel; } 
  void setSDiel( double sdiel ) { m_sdiel = sdiel; resetKappaFact2( ); } 
  void setTemp( double temp ) { m_temp = temp; resetKappaFact2( ); }
  void setRand( int rand ) { m_srand = rand; }
//  void setDV( double dVolt ) { m_fExt[0] = dVolt; }
//  void setMThick( double mThick ) { m_fExt[1] = mThick; } 
//  void setFExt( double dVolt, double mthick )
//                { setDV(dVolt); setMThick(mthick); }
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
//  void setTypeNImat( int typeCount, string imat ) 
//              { m_molfnames[typeCount][1] = imat; }
//  void setTypeNSpolDir( int typeCount, string poldir ) 
//              { m_molfnames[typeCount][2] = poldir; }
//  void setTypeNExp( int typeCount, string ext ) 
//              { m_molfnames[typeCount][3] = ext; }
  void setTypeNXYZ( int typeCount, string xyz ) 
              { m_molfnames[typeCount][1] = xyz; }
  
//  void setScale( double scal ) { m_scale = scal; }
//  void setMolDist( double dd ) { m_molD = dd; }
//  void setLayer( int layer ) { m_layer = layer; }
	
	string getRunType()              { return m_runSpecs[0]; }
//	string getSimType()              { return m_runSpecs[1]; }
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

//  double getDist()                 { return m_distance; }

//  double getDVolt()                { return m_fExt[0]; }
//  double getMThick()               { return m_fExt[1]; }
  int getNTraj()                   { return m_ntraj; }

//  double getCutoffInter() { return m_cutInter ; }
//  double getCutoffInteract() { return m_cutInteract ; }
//  double getCutoffIntra() { return m_cutIntra ; }

  double getDtr( int n )           { return m_typeDiff[n][0]; }
  double getDrot( int n )          { return m_typeDiff[n][1]; }

	int getTypeNCount(int type)      { return m_nTypenCount[type]; }
	string getTypeNDef(int type)     { return m_typeDef[type]; }
	string getTypeNPQR(int type)     { return m_molfnames[type][0]; }
//	string getTypeNMatDir(int type)  { return m_molfnames[type][1]; }
//	string getTypeNPolDir(int type)  { return m_molfnames[type][2]; }
//	string getTypeNPolExp(int type)  { return m_molfnames[type][3]; }
	string getTypeNXYZ(int type)     { return m_molfnames[type][1]; }

  double getKappa()                { return m_kappa; }
  double getIKbT()                 { return m_IKbT; }
  double getFACT2()                { return m_fact2; }
//  double getScale()                { return m_scale; }
//  double getMolD()                 { return m_molD; }
//  double getLayer()                { return m_layer; }

private:
  Constants consts_;
  int m_ompThreads;
  double m_saltConc;
  int m_nType;  		// Number of different molecule types
  int m_PBCs;			// PBC in 0, pseudo-2, or 3 dimensions
  double m_blen; 		// boxlength for PBC
  int m_maxtime;

//  double m_distance;  // distance between bodies

//  double m_cutInter , m_cutInteract ;
//  double m_cutIntra ; // cut-offs for inter- intra-molecular interact in Ang

  int m_ntraj;

  double m_idiel, m_sdiel;// dielectric constant win molecule and of solvent
  double m_temp; 
  int m_srand;			// random seed

//  double m_fExt[2];	//External potential drop params, 0=d volt, 1=memb thick
  vector<int> m_nTypenCount; // Array for each of mol types, how many mols
  vector<vector<double> > m_typeDiff; // Dtr,Drot  each type, size [Ntype][2]
  vector<string> m_typeDef; 		// For each type, type is stat, rot or move
  string m_runSpecs[2];	//include run type [0],
                        // pbsam/bd and the runname [1]
  vector<vector<string> > m_molfnames;  // file names 

  double m_kappa, m_IKbT, m_fact2;
//  m_scale, m_molD;
//  int m_layer; // The number of layers for an inf grid
}; // end CSetup

#endif

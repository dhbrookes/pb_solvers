#include "setup.h"

/*#########################################################*/
/*#########################################################*/
// CSetup constructor
/*#########################################################*/
/*#########################################################*/
CSetup::CSetup(Constants consts) :
  m_ompThreads( 1 ),
  m_saltConc( 0.01 ),
  m_nType( 2 ),
  m_PBCs( 0 ),
  m_blen( 1.4e8 ),
  m_maxtime( 1000000 ),
  m_ntraj ( 5 ),
  m_idiel( 4.0 ),
  m_sdiel( 78.0 ),
  m_temp( 298.0 ),
  m_srand( (unsigned)time(NULL) ),
//  m_cutInter( 10 ) ,
//  m_cutInteract( 100 ) ,
//  m_cutIntra( 30 ),
//  m_scale( 1.0 ),
//  m_molD( 100.0 ),
//  m_layer( 1 ),
  consts_(consts)
{
//  m_fExt[0] = 0.0; m_fExt[1] =  1.0;
  m_nTypenCount.resize(m_nType); m_typeDef.resize(m_nType);
  m_nTypenCount[0] = 1; m_nTypenCount[1] = 1;

  m_typeDiff.resize(m_nType);
  for (int i = 0; i<m_nType; i++)
  {
    m_typeDiff[i].resize(2);
    m_typeDiff[i][0] = 0.0; m_typeDiff[i][1] = 0.0;
  }

  m_typeDef[0] = "stat"; m_typeDef[1] = "stat";
  m_runSpecs[0] = "pot"; m_runSpecs[1] = "nafion"; m_runSpecs[2] = "test";

  // Initializing file locs to defaults
  // pqr fname, imat path, spol path, spol name
  const char * molfns[][5] = {{"../Config/test1.pqr", "../Imat/test1/",
                  "../Selfpol/test1", "test1_p30.0", "../Config/test1.xyz"},
                  {"../Config/test2.pqr", "../Imat/test2/",
                  "../Selfpol/test2", "test2_p30.0", "../Config/test2.xyz"}};

  m_molfnames.resize(m_nType);
  for (int i=0; i<m_nType; i++)
  {
    m_molfnames[i].resize(5);
    for (int j=0; j<m_molfnames[i].size();j++)
      m_molfnames[i][j] = molfns[i][j];
  }

  resetKappaFact2();
} // end CSetup constructor


/*#########################################################*/
/*#########################################################*/
// resizeVecs
// A function to resize necessary vecs in CSetup class
/*#########################################################*/
/*#########################################################*/
void CSetup::resizeVecs()
{
  m_nTypenCount.resize(m_nType);
  m_typeDef.resize(m_nType);

  m_typeDiff.resize(m_nType);
  for (int i = 0; i<m_nType; i++) m_typeDiff[i].resize(2);

  m_molfnames.resize(m_nType);
  for(int i = 0; i < m_nType; i++)
    m_molfnames[i].resize(5);
} // end resizeVecs

/*#########################################################*/
/*#########################################################*/
// printSetupClass
// A function to print variables of the setup class
/*#########################################################*/
/*#########################################################*/
void CSetup::printSetupClass( )
{
  cout << "This are my setup parameters: " << endl;
  cout << "kappa: " << m_kappa << " temp: " << m_temp <<
     " IKbT: " << m_IKbT << " Fact2: " << m_fact2 <<
     " idiel: " << m_idiel << " sdiel: " << m_sdiel << endl;
  cout << "N OMP: " << m_ompThreads << " Salt con: " << m_saltConc <<
     " N mol types: " << m_nType<< " Max time: " << m_maxtime << endl;
  cout << " N Trajs: " << m_ntraj<< " PBC type: "  << m_PBCs <<
     " Box len: " << m_blen << " Random Seed: " <<  m_srand<< endl;
  cout << " Run specs, sys: " << m_runSpecs[0] << "; Type: " <<
     m_runSpecs[1] << "; Runname: " << m_runSpecs[2] << endl;


  for (int i = 0; i<m_nType; i++)
  {
    cout << "MolType: " << i << " count: " << m_nTypenCount[i];
    cout << " Def: " << m_typeDef[i] << ": Dtr: ";
    cout << m_typeDiff[i][0] << ": Drot: " << m_typeDiff[i][1] << endl;
  }

  cout << endl;
} // end printSetupClass

/*#########################################################*/
/*#########################################################*/
// resetKappa
// A function to reset kappa and other random system parameters
/*#########################################################*/
/*#########################################################*/
void CSetup::resetKappaFact2( )
{
  m_kappa = consts_.ANGSTROM * sqrt( (2*m_saltConc*
                   consts_.AVOGADRO_NUM/consts_.LITRE*consts_.E2 )
                   / ( m_sdiel* consts_.PERMITTIVITY_VAC * consts_.KB * m_temp  ) );

  m_IKbT = 1.0/( consts_.KB*m_temp ); // [1/J]
  m_fact2 = consts_.convert_int_to_kcal_mol(m_IKbT);
} //end resetKappaFact2


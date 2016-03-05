#include "setup.h"

/*#########################################################*/
/*#########################################################*/
// CSetup constructor
/*#########################################################*/
/*#########################################################*/
CSetup::CSetup() :
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
  m_srand( (unsigned)time(NULL) )
{
  m_nTypenCount.resize(m_nType); m_typeDef.resize(m_nType);
  m_nTypenCount[0] = 1; m_nTypenCount[1] = 1;

  m_typeDiff.resize(m_nType);
  for (int i = 0; i<m_nType; i++)
  {
    m_typeDiff[i].resize(2);
    m_typeDiff[i][0] = 0.0; m_typeDiff[i][1] = 0.0;
  }

  m_typeDef[0] = "stat"; m_typeDef[1] = "stat";
  m_runSpecs[0] = "pot"; m_runSpecs[2] = "test";

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
     " IKbT: " << m_IKbT <<
  " idiel: " << m_idiel << " sdiel: " << m_sdiel << endl;
  cout << "N OMP: " << m_ompThreads << " Salt con: " << m_saltConc <<
     " N mol types: " << m_nType<< " Max time: " << m_maxtime << endl;
  cout << " N Trajs: " << m_ntraj<< " PBC type: "  << m_PBCs <<
     " Box len: " << m_blen << " Random Seed: " <<  m_srand<< endl;
//  cout << " F Ext, DV: " << m_fExt[0] << " M thick: " << m_fExt[1] << endl;
  cout << " Run specs, sys: " << m_runSpecs[0] << "; Type: " <<
     m_runSpecs[1] << "; Runname: " << m_runSpecs[2] << endl;
//  cout << " Cutoffs : Interact = " << m_cutInteract << " ; inter = " <<
//     m_cutInter << " ; intra = " << m_cutIntra << endl ; 

  for (int i = 0; i<m_nType; i++)
  {
    cout << "MolType: " << i << " count: " << m_nTypenCount[i];
    cout << " Def: " << m_typeDef[i] << ": Dtr: ";
    cout << m_typeDiff[i][0] << ": Drot: " << m_typeDiff[i][1] << endl;
  }

  cout << endl;
} // end printSetupClass

Constants CSetup::setup_constants()
{
  Constants cons = Constants();
  cons.set_dielectric_prot(m_idiel);
  cons.set_dielectric_water(m_sdiel);
  cons.set_salt_concentration(m_saltConc);
  cons.set_temp(m_temp);
  cons.update_all();
  return cons;
}

System CSetup::setup_system(Constants consts)
{
  vector<Molecule> mols;
  
  int i, j;
  string pqrpath;
  for (i = 0; i < m_nType; i++)
  {
    for (j = 0; j < m_nTypenCount[i]; j++)
    {
      pqrpath = m_molfnames[i][j];
      PQRFile pqrobj (pqrpath);
      Molecule mol (m_typeDef[i], pqrobj.get_charges(), pqrobj.get_pts(),
               pqrobj.get_radii(), getDrot(i), getDtr(i));
      mols.push_back(mol);
    }
  }
  System sys = System(consts, mols);
  return sys;
}


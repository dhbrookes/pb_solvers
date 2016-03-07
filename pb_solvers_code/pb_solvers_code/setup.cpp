
#include "setup.h"

///*#########################################################*/
///*#########################################################*/
//// CSetup constructor
///*#########################################################*/
///*#########################################################*/
//CSetup::CSetup() :
//  m_ompThreads( 1 ),
//  m_saltConc( 0.01 ),
//  m_nType( 2 ),
//  m_PBCs( 0 ),
//  m_blen( 1.4e8 ),
//  m_maxtime( 1000000 ),
//  m_ntraj ( 5 ),
//  m_idiel( 4.0 ),
//  m_sdiel( 78.0 ),
//  m_temp( 298.0 ),
//  m_srand( (unsigned)time(NULL) )
//{
//  m_nTypenCount.resize(m_nType); m_typeDef.resize(m_nType);
//  m_nTypenCount[0] = 1; m_nTypenCount[1] = 1;
//
//  m_typeDiff.resize(m_nType);
//  for (int i = 0; i<m_nType; i++)
//  {
//    m_typeDiff[i].resize(2);
//    m_typeDiff[i][0] = 0.0; m_typeDiff[i][1] = 0.0;
//  }
//
//  m_typeDef[0] = "stat"; m_typeDef[1] = "stat";
//  m_runSpecs[0] = "pot"; m_runSpecs[2] = "test";
//
//  // Initializing file locs to defaults
//  // pqr fname, imat path, spol path, spol name
//  const char * molfns[][5] = {{"../Config/test1.pqr", "../Imat/test1/",
//                  "../Selfpol/test1", "test1_p30.0", "../Config/test1.xyz"},
//                  {"../Config/test2.pqr", "../Imat/test2/",
//                  "../Selfpol/test2", "test2_p30.0", "../Config/test2.xyz"}};
//
//  m_molfnames.resize(m_nType);
//  for (int i=0; i<m_nType; i++)
//  {
//    m_molfnames[i].resize(5);
//    for (int j=0; j<m_molfnames[i].size();j++)
//      m_molfnames[i][j] = molfns[i][j];
//  }
//} // end CSetup constructor
//
//
///*#########################################################*/
///*#########################################################*/
//// resizeVecs
//// A function to resize necessary vecs in CSetup class
///*#########################################################*/
///*#########################################################*/
//void CSetup::resizeVecs()
//{
//  m_nTypenCount.resize(m_nType);
//  m_typeDef.resize(m_nType);
//
//  m_typeDiff.resize(m_nType);
//  for (int i = 0; i<m_nType; i++) m_typeDiff[i].resize(2);
//
//  m_molfnames.resize(m_nType);
//  for(int i = 0; i < m_nType; i++)
//    m_molfnames[i].resize(5);
//} // end resizeVecs
//
///*#########################################################*/
///*#########################################################*/
//// printSetupClass
//// A function to print variables of the setup class
///*#########################################################*/
///*#########################################################*/
//void CSetup::printSetupClass( )
//{
//  cout << "This are my setup parameters: " << endl;
//  cout << "kappa: " << m_kappa << " temp: " << m_temp <<
//     " IKbT: " << m_IKbT <<
//  " idiel: " << m_idiel << " sdiel: " << m_sdiel << endl;
//  cout << "N OMP: " << m_ompThreads << " Salt con: " << m_saltConc <<
//     " N mol types: " << m_nType<< " Max time: " << m_maxtime << endl;
//  cout << " N Trajs: " << m_ntraj<< " PBC type: "  << m_PBCs <<
//     " Box len: " << m_blen << " Random Seed: " <<  m_srand<< endl;
////  cout << " F Ext, DV: " << m_fExt[0] << " M thick: " << m_fExt[1] << endl;
//  cout << " Run specs, sys: " << m_runSpecs[0] << "; Type: " <<
//     m_runSpecs[1] << "; Runname: " << m_runSpecs[2] << endl;
////  cout << " Cutoffs : Interact = " << m_cutInteract << " ; inter = " <<
////     m_cutInter << " ; intra = " << m_cutIntra << endl ; 
//
//  for (int i = 0; i<m_nType; i++)
//  {
//    cout << "MolType: " << i << " count: " << m_nTypenCount[i];
//    cout << " Def: " << m_typeDef[i] << ": Dtr: ";
//    cout << m_typeDiff[i][0] << ": Drot: " << m_typeDiff[i][1] << endl;
//  }
//
//  cout << endl;
//} // end printSetupClass
//
//Constants CSetup::setup_constants()
//{
//  Constants cons = Constants();
//  cons.set_dielectric_prot(m_idiel);
//  cons.set_dielectric_water(m_sdiel);
//  cons.set_salt_concentration(m_saltConc);
//  cons.set_temp(m_temp);
//  cons.update_all();
//  return cons;
//}
//
//System CSetup::setup_system(Constants consts)
//{
//  vector<Molecule> mols;
//  
//  int i, j;
//  string pqrpath;
//  for (i = 0; i < m_nType; i++)
//  {
//    for (j = 0; j < m_nTypenCount[i]; j++)
//    {
//      pqrpath = m_molfnames[i][j];
//      PQRFile pqrobj (pqrpath);
//      Molecule mol (m_typeDef[i], pqrobj.get_charges(), pqrobj.get_pts(),
//               pqrobj.get_radii(), getDrot(i), getDtr(i));
//      mols.push_back(mol);
//    }
//  }
//  System sys = System(consts, mols);
//  return sys;
//}


Setup::Setup(string infile)
:ompThreads_( 1 ),
saltConc_( 0.01 ),
nType_( 2 ),
PBCs_( 0 ),
blen_( 1.4e8 ),
maxtime_( 1000000 ),
ntraj_( 5 ),
idiel_( 4.0 ),
sdiel_( 78.0 ),
temp_( 298.0 ),
srand_( (unsigned)time(NULL) ),
nTypenCount_(2),
typeDef_(2),
typeDiff_(2),
molfnames_(2),
runSpecs_(2)
{
  nTypenCount_[0] = 1;
  nTypenCount_[1] = 1;
  
  for (int i = 0; i<nType_; i++)
  {
    typeDiff_[i] = vector<double> (2);
    typeDiff_[i][0] = 0.0;
    typeDiff_[i][1] = 0.0;
  }
  
  typeDef_[0] = "stat";
  typeDef_[1] = "stat";
  runSpecs_[0] = "pot";
  runSpecs_[1] = "test";
  
  // Initializing file locs to defaults
  // pqr fname, imat path, spol path, spol name
  vector<vector<string> > molfn = {{"../Config/test1.pqr", "../Imat/test1/",
    "../Selfpol/test1", "test1_p30.0", "../Config/test1.xyz"},
    {"../Config/test2.pqr", "../Imat/test2/",
      "../Selfpol/test2", "test2_p30.0", "../Config/test2.xyz"}};
  
  for (int i=0; i<nType_; i++)
  {
    molfnames_[i] = vector<string> ( 5);
    for (int j=0; j<molfnames_[i].size();j++)
    {
      molfnames_[i][j] = molfn[i][j];
    }
  }
  
}


void Setup::read_infile(string fname)
{
  cout << "Reading Input file " << fname << endl ;
  ifstream fin(fname);
  if (!fin.is_open())
  {
    cout << "Could not open file " << fname << endl;
    exit(0);
  }
  
  string inputLine;
  vector<vector <string> > keywordLines;
  getline(fin,inputLine);
  
  while (!fin.eof())
  {
    findLines(inputLine);
    getline(fin, inputLine);
  }
}

vector<string> Setup::split(string str, char delim)
{
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  
  while(getline(ss, tok, delim))
  {
    internal.push_back(tok);
  }
  return internal;
}


void Setup::findLines(string fline)
{
    if (!fline.empty())
    {
      findKeyword( split(fline, ' '));
    }
}

void Setup::findKeyword(vector<string> fline)
{
  string keyword = fline[0];
  if (keyword == "runname")
  {
    cout << "Runname command found" << endl;
    setRunName( fline[1] );
  } else if (keyword == "runtype")
  {
    cout << "Runtype command found" << endl;
    setRunType( fline[1] );
    if (fline[1] == "sim")
    {
      if (fline.size() > 2)
      {
        setNTraj( atoi( fline[2].c_str() ) );
      }
      if (fline.size() > 3)
      {
        setMaxTime( atoi( fline[3].c_str() ));
      }
    }
  } else if (keyword == "omp")
  {
    cout << "OMP command found" << endl;
    setOMP(atoi(fline[1].c_str()));
  } else if (keyword == "pbc")
  {
    cout << "PBC command found" << endl;
    setPBCT( atoi(fline[1].c_str()) );
    if ( getPBCs() > 0 )
      setBoxl(atof(fline[2].c_str()));
  } else if (keyword == "salt")
  {
    cout << "Salt command found" << endl;
    setSaltCon( atof(fline[1].c_str()) );
  } else if (keyword == "temp")
  {
    cout << "Temperature command found" << endl;
    setTemp( atof(fline[1].c_str()) );
  } else if (keyword == "idiel")
  {
    cout << "Interior dielectric command found" << endl;
    setIDiel( atof(fline[1].c_str()) );
  } else if (keyword == "sdiel")
  {
    cout << "Solvent dielectric command found" << endl;
    setSDiel( atof(fline[1].c_str()) );
  } else if (keyword == "attypes")
  {
    cout << "Atom Types command found" << endl;
    setNType( atoi(fline[1].c_str()) );
    resizeVecs();
    cout << "done with attypes" << endl;
  } else if (keyword == "type")
  {
    cout << "Type def command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getType()-1)
      return;
    if (fline.size() > 2)
      setTypeNCount( typeNo, atoi(fline[2].c_str()) );
    if (fline.size() > 3)
    {
      setTypeNDef( typeNo, fline[3].c_str() );
      if (getTypeNDef(typeNo) == "move")
        setTypeNDtr( typeNo, atof(fline[4].c_str()));
    }
  } else if (keyword == "pqr")
  {
    cout << "PQR command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getType()-1)
      return;
    setTypeNPQR( typeNo, fline[2].c_str() );
  } else if (keyword == "xyz")
  {
    cout << "XYZ command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getType()-1)
      return;
    setTypeNXYZ( typeNo, fline[2] );
    
  } else if (keyword == "random")
  {
    cout << "RNG Seed command found" << endl;
    setRand( atoi(fline[1].c_str()) );
  } else
    cout << "Keyword not found, read in as " << fline[0] << endl;
}

void Setup::resizeVecs()
{
  nTypenCount_.resize(nType_);
  typeDef_.resize(nType_);

  typeDiff_.resize(nType_);
  for (int i = 0; i<nType_; i++)
  {
    typeDiff_[i].resize(2);
  }

  molfnames_.resize(nType_);
  for(int i = 0; i < nType_; i++)
  {
    molfnames_[i].resize(5);
  }
} // end resizeVecs


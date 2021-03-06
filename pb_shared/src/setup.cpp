//
//  setup.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 3/9/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "setup.h"

const double Setup::MAX_DIST = 1.4e18;

Setup::Setup(string infile)
:
ompThreads_( 1 ),
saltConc_( 0.01 ),
nType_( 2 ),
PBCs_( 0 ),
blen_( MAX_DIST ),
maxtime_( 1000000 ),
ntraj_( 5 ),
gridPts_( 30 ),
gridCt_(0),
idiel_( 4.0 ),
sdiel_( 78.0 ),
temp_( 298.0 ),
npoles_( 5 ),
srand_( (unsigned)time(NULL) ),
nTypenCount_(2),
typeDef_(2),
typeDiff_(2),
pqr_names_(2),
surfNames_(2),
imatNames_(2),
expNames_(2),
xyz_names_(2),
isTransRot_(2),
runSpecs_(2),
mbdfile_loc_(2),
termvals_(2),
termtype_(2),
andCombine_(false),
sphBeta_(2.0),
tolSP_(1.0),
maxTrials_(40),
nTrials_(9)
{
  nTypenCount_[0] = 1;
  nTypenCount_[1] = 1;
  
  for (int i = 0; i<nType_; i++)
  {
    typeDiff_[i] = vector<double> (2);
    typeDiff_[i][0] = 0.0;
    typeDiff_[i][1] = 0.0;
  }
  
  typeDef_[0]  = "stat";
  typeDef_[1]  = "stat";
  runSpecs_[0] = "electrostatics";
  runSpecs_[1] = "test";
  
  potOutfnames_.resize(3);
  potOutfnames_[0] = "";
  potOutfnames_[1] = "";
  potOutfnames_[2] = "";
  
  mbdfile_loc_[0] = "";
  mbdfile_loc_[1] = "";
  
  imatNames_[0] = "";
  imatNames_[1] = "";
  
  expNames_[0] = "";
  expNames_[1] = "";
  
  units_ = "internal";
  
  // Initializing file locs to defaults
  // pqr fname, imat path, spol path, spol name
  vector<vector<string> > molfn = {{"../Config/test1.pqr",
    "../Config/test1.xyz"}, {"../Config/test2.pqr", "../Config/test2.xyz"}};
  
  for (int i=0; i<nType_; i++)
  {
    pqr_names_[i] = molfn[i][0];
    xyz_names_[i].resize(1);
    isTransRot_[i].resize(1);
    xyz_names_[i][0] = molfn[i][1];
    isTransRot_[i][0] = false;
  }
  
  confiles_.resize(0);
  
  read_infile(infile);
}

// APBS
Setup::Setup(double temp, double salt_conc, double int_diel, double solv_diel,
             int nmol, string runtype, string runname, bool randorient, 
             double boxl, int pbc_type, int gridpts, string map3d, int g2dct, 
             vector<string> grid2Dfn, vector <string> grid2Dax, 
             vector <double> grid2Dloc, string dxnam, int ntraj, bool termcomb,
             vector<string> difftype, vector<vector<double> > diffcon,
             vector<string> termcond, vector<double> termval, 
             vector<vector <int > > termnu, vector<string> confil,
             vector<double> conpad, vector<vector <string> > xyzfil, string unt)
:
ompThreads_( 1 ),
saltConc_( salt_conc ), //
nType_( nmol ),  //
PBCs_( pbc_type ),  //
blen_( boxl ),   //
maxtime_( 1000000 ),
ntraj_( ntraj ), //
gridPts_( gridpts ), //
gridCt_(g2dct), //
idiel_( int_diel ),  //
sdiel_( solv_diel ), //
temp_( temp ),       //
srand_( (unsigned)time(NULL) ),
nTypenCount_(nmol), //
typeDef_(nmol),
typeDiff_(nmol),
pqr_names_(nmol),
xyz_names_(nmol),
imatNames_(nmol),
surfNames_(nmol),
expNames_(nmol),
isTransRot_(nmol),
runSpecs_(2),  // 
mbdfile_loc_(2),
numTerm_((int)termcond.size()), //
termvals_((int)termcond.size()), //
termmols_((int)termcond.size()), //
termtype_((int)termcond.size()), //
confiles_((int)confil.size()),  //
conpads_((int)confil.size()),  //
andCombine_(termcomb), //
orientRand_(randorient) //
{
  runSpecs_[0] = runtype; //
  runSpecs_[1] = runname; //
  units_ = unt;

  for (int i = 0; i<nType_; i++) nTypenCount_[i] = 1; //

  // Dynamics part
  confiles_.resize(0);
  for (int i = 0; i<nType_; i++)
  {
    typeDef_[i] = difftype[i];
    typeDiff_[i] = vector<double> (2);
    typeDiff_[i][0] = diffcon[i][0];
    typeDiff_[i][1] = diffcon[i][1];
  }

  for (int i = 0; i < numTerm_; i++)
  {
    termtype_[i] = termcond[i];
    termmols_[i] = vector<int> (1);
    termmols_[i][0] = termnu[i][0];
    termvals_[i] = termval[i];
  }

  for (int i = 0; i < confiles_.size(); i++)
  {
    confiles_[i] = confil[i];
    conpads_[i] = conpad[i];
  }

  // Electrostatics part
  potOutfnames_.resize(2+g2dct); //
  axis_.resize(g2dct);  //
  axLoc_.resize(g2dct); //
  potOutfnames_[0] = dxnam;  //
  potOutfnames_[1] = map3d;  //
  
  for (int i = 0; i<grid2Dax.size(); i++) 
  {
    potOutfnames_[2+i] = grid2Dfn[i]; //
    axis_[i] = grid2Dax[i];  //
    axLoc_[i] = grid2Dloc[i];
  }
  
  // Mutibody expansion part
  mbdfile_loc_[0] = "";
  mbdfile_loc_[1] = "";

  // Molecule part
  for (int i=0; i<nType_; i++)
  {
    pqr_names_[i] = "mol" + to_string(i) + ".pqr";
    surfNames_[i] = "";
    imatNames_[i] = "";
    expNames_[i]  = "";
    xyz_names_[i].resize(ntraj);
    isTransRot_[i].resize(1);
    isTransRot_[i][0] = false;

    for (int j=0; j<ntraj; j++) xyz_names_[i][j] = xyzfil[i][j];
  }  
}

void Setup::apbs_pbsam_set(vector<string> surffil, vector<string> imatfil,
                           vector<string> expfil)
{
  int i;
//if ( nType_ > surffil.size() ) 
//  cout << "Missing surf file, assuming its the end ones" << endl;
  for (i=0; i<surffil.size(); i++)
    surfNames_[i] = surffil[i];

//if ( nType_ > imatfil.size() ) 
//  cout << "Missing imat file, assuming its the end ones" << endl;
  for (i=0; i<imatfil.size(); i++)
    imatNames_[i] = imatfil[i];

//if ( nType_ > expfil.size() ) 
//  cout << "Missing exp file, assuming its the end ones" << endl;
  for (i=0; i<expfil.size(); i++)
    expNames_[i] = expfil[i];
}

void Setup::read_infile(string fname)
{
  cout << "Reading Input file " << fname << endl ;
  ifstream fin(fname);
  if (!fin.is_open()) throw CouldNotReadException(fname);
  
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
    if (fline[1] == "dynamics")
    {
      if (fline.size() > 2)
      {
        setNTraj( atoi( fline[2].c_str() ) );
      }
      if (fline.size() > 3)
      {
        setMaxTime( atoi( fline[3].c_str() ));
      }
    } else if (fline[1] == "electrostatics")
    {
      if (fline.size() > 2)
      {
        setGridPts( atoi( fline[2].c_str() ));
      }
    }
  } else if (keyword == "dx")
  {
    cout << "DX command found" << endl;
    setDXoutName( fline[1].c_str());
  } else if (keyword == "3dmap")
  {
    cout << "3D map command found" << endl;
    set_3dmap_name( fline[1].c_str());
  } else if (keyword == "gridct")
  {
    cout << "Grid count command found" << endl;
    setGridCt( atoi(fline[1].c_str()));
    axis_.resize( gridCt_); axLoc_.resize(gridCt_);
    potOutfnames_.resize(gridCt_+2);
  } else if (keyword == "grid2D")
  {
    cout << "Grid command found" << endl;
    setGridOutName( atoi(fline[1].c_str()), fline[2].c_str());
    setGridAx( atoi(fline[1].c_str()), fline[3].c_str());
    setGridAxLoc( atoi(fline[1].c_str()), atof(fline[4].c_str()));
  } else if (keyword == "3bdloc")
  {
    cout << "3BD loc command found" << endl;
    set3BDLoc(fline[1].c_str());
  } else if (keyword == "2bdloc")
  {
    cout << "2BD loc command found" << endl;
    set2BDLoc(fline[1].c_str());
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
  } else if (keyword == "poles")
  {
    cout << "Number of poles command found" << endl;
    setNPoles( atoi(fline[1].c_str()) );
  }else if (keyword == "termct")
  {
    cout << "Termination count command found" << endl;
    set_numterms(atoi(fline[1].c_str()));
    resize_termcond(atoi(fline[1].c_str()));
  } else if (keyword == "termcombine")
  {
    cout << "Termination combine command found" << endl;
    set_term_combine(fline[1]);
  } else if (keyword == "term")
  {
    cout << "Termination condition command found" << endl;
    int idx = atoi(fline[1].c_str()) - 1;
    if (idx > get_numterms()-1)
    {
      cout << "WARNING: trying to add more term types than specified" << endl;
      return;
    }
    string type = fline[2];
    double val = atof(fline[3].c_str());
    vector<int> mol_idx(2);
    if (type == "contact")
    {
      string confile = fline[3];
      double pad  = atof(fline[4].c_str());
      // placeholders:
      mol_idx = {-1};
      val = -1;
      confiles_.push_back(fline[3]);
      cout << "Contact size " << confiles_.size() << endl;
      conpads_.push_back(pad);
      cout << "This is my first contact file " << confiles_[0] << endl;
      cout << "This is my first contact file " << fline[3] << endl;
    }
    else
    {
      mol_idx[0] = atoi(fline[4].c_str()) - 1;
    }
    add_termcond(idx, type, mol_idx, val);
    
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
    if (typeNo > getNType()-1)
    {
      cout << "WARNING: trying to add more mol types than specified" << endl;
      return;
    }
    if (fline.size() > 2)
      setTypeNCount( typeNo, atoi(fline[2].c_str()) );
    if (fline.size() > 3)
    {
      setTypeNDef( typeNo, fline[3].c_str() );
      if (getTypeNDef(typeNo) == "move")
      {
        setTypeNDtr( typeNo, atof(fline[4].c_str()));
        setTypeNDrot( typeNo, atof(fline[5].c_str()));
      } else if (getTypeNDef(typeNo) == "rot")
      {
        setTypeNDtr( typeNo, 0.0);
        setTypeNDrot( typeNo, atof(fline[4].c_str()));
      } else
      {
        setTypeNDtr( typeNo, 0.0);
        setTypeNDrot( typeNo, 0.0);
      }
    }
  } else if (keyword == "pqr")
  {
    cout << "PQR command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getNType()-1)
      return;
    setTypeNPQR( typeNo, fline[2].c_str() );
  } else if (keyword == "xyz")
  {
    string xyz;
    int traj, typeNo = atoi(fline[1].c_str())-1;
    
    if ( fline.size() == 4 )
    {
      traj = atoi(fline[2].c_str())-1;
      xyz = fline[3];
    } else
    {
      traj = 0;
      xyz = fline[2];
    }
    if (typeNo > getNType()-1)
      return;
    setTypeNXYZ( typeNo, traj, xyz );
    cout << "XYZ command found " << xyz << endl;
    setTypeNisTransRot(typeNo, traj, false);
  } else if (keyword == "transrot")
  {
    string transrot;
    int traj, typeNo = atoi(fline[1].c_str())-1;
    cout << "transrot command found" << endl;
    
    if ( fline.size() == 4 )
    {
      traj = atoi(fline[2].c_str())-1;
      transrot = fline[3];
    } else
    {
      traj = 0;
      transrot = fline[2];
    }
    if (typeNo > getNType()-1)
      return;
    setTypeNXYZ( typeNo, traj, transrot );
    setTypeNisTransRot(typeNo, traj, true);
  } else if (keyword == "surf")
  {
    cout << "surf command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getNType()-1)
      return;
    setTypeNSurf( typeNo, fline[2].c_str() );
  } else if (keyword == "imat")
  {
    cout << "IMAT prefix command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getNType()-1)
      return;
    setTypeNImat( typeNo, fline[2].c_str() );
  } else if (keyword == "exp")
  {
    cout << "Expansion prefix command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > getNType()-1)
      return;
    setTypeNExp( typeNo, fline[2].c_str() );
  } else if (keyword == "randorient")
  {
    cout << "Random orientation command found" << endl;
    setRandOrient();
  } else if (keyword == "random")
  {
    cout << "RNG Seed command found" << endl;
    setRand( atoi(fline[1].c_str()) );
  } else if (keyword == "units")
  {
    cout << "Units command found" << endl;
    setUnits( fline[1].c_str() );
  } else if (keyword == "tolsp")
  {
    cout << "tolsp command found" << endl;
    set_tol_sp(atof(fline[1].c_str()));
  } else if (keyword == "ntrials")
  {
    cout << "ntrials command found" << endl;
    set_n_trials(atoi(fline[1].c_str()));
  } else if (keyword == "maxtrials")
  {
    cout << "maxtrials command found" << endl;
    set_max_trials(atoi(fline[1].c_str()));
  } else if (keyword == "sphbeta")
  {
    cout << "sphbeta command found" << endl;
    set_sph_beta(atof(fline[1].c_str()));
  } else
    cout << "Keyword not found, read in as " << fline[0] << endl;
}

void Setup::resizeVecs()
{
  int ntraj= (getRunType() == "dynamics") ? getNTraj() : 1;
  nTypenCount_.resize(nType_);
  typeDef_.resize(nType_);

  typeDiff_.resize(nType_);
  for (int i = 0; i<nType_; i++)
  {
    typeDiff_[i].resize(2);
  }

  pqr_names_.resize(nType_);
  xyz_names_.resize(nType_);
  surfNames_.resize(nType_);
  imatNames_.resize(nType_);
  expNames_.resize(nType_);
  isTransRot_.resize(nType_);
  for(int i = 0; i < nType_; i++)
  {
    xyz_names_[i].resize(ntraj);
    surfNames_[i] = "";
    imatNames_[i] = "";
    expNames_[i] = "";
    isTransRot_[i].resize(nType_);
  }

} // end resizeVecs

void Setup::check_inputs()
{
  vector<string> problems;
  if (typeDiff_.size() < nType_ )
  {
    problems.push_back("Number of molecular input type parameters is less \
                       than the specified number of MoleculeSAMAM types");
  }
  if (typeDef_.size() < nType_)
  {
    problems.push_back("Number of movement types is less \
                        than the specified number of MoleculeSAMAM types");
  }
  if (pqr_names_.size() < nType_)
  {
    problems.push_back("Number of provided PQR files is less than the\
                       specified number of molecular types");
  }
  if (xyz_names_.size() < nType_)
  {
    problems.push_back("Number of provided XYZ files is less than the\
                       specified number of molecular types");
    for (int i = 0; i < nType_; i++)
      if (xyz_names_[i].size() < ntraj_)
      {
        problems.push_back("Number of provided XYZ files is less than the\
                           specified number of trajectories");
      }
  }
  if (runSpecs_[0] == "electrostatics")
  {
    if (axis_.size() < gridCt_ || axLoc_.size() < gridCt_)
    {
      problems.push_back("Number grids provided is less than specified \
                          grid count");
    }
  }
  if (runSpecs_[0] == "dynamics")
  {
    if (termtype_.size() < numTerm_ || termtype_.size() < numTerm_)
    {
      problems.push_back("Number of termination conditions provided is less \
                              than specified termination count");
    }
  }
  if (problems.size() > 0) throw BadInputException(problems);
}


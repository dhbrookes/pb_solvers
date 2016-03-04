#include "readinput.h"

/*!#########################################################*/
/*#########################################################*/
// split
// Split a input line and return a vector
/*#########################################################*/
/*#########################################################*/
vector<string> split(string str, char delimiter)
{
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;

  while(getline(ss, tok, delimiter))
    internal.push_back(tok);

  return internal;
} // end split


/*!#########################################################*/
/*#########################################################*/
// findLines
// Read an input line and break up if not blank
/*#########################################################*/
/*#########################################################*/
void findLines(string fline, CSetup & inputRead)
{
  if (!fline.empty())
    findKeyword( split(fline, ' '), inputRead);

} // end findLines

/*!#########################################################*/
/*#########################################################*/
// findKeyword
// Read an input line and determine relevant parameters
/*#########################################################*/
/*#########################################################*/
void findKeyword(vector<string> fline, CSetup & inputRead)
{
  string keyword = fline[0];

//  if (keyword == "system")
//  {
//    cout << "System command found" << endl;
//    if(fline[1]=="nafion" || fline[1]=="multi" || fline[1] == "nam" )
//    {
//      inputRead.setSimType( fline[1] );
//    }
//    else
//       cout << "Option not available, using default" << endl;

  if (keyword == "runname")
  {
    cout << "Runname command found" << endl;
    inputRead.setRunName( fline[1] );
  } else if (keyword == "runtype")
  {
    cout << "Runtype command found" << endl;
    inputRead.setRunType( fline[1] );
    if (fline[1] == "sim")
    {
      if (fline.size() > 2)
        inputRead.setNTraj( atoi( fline[2].c_str() ) );
      if (fline.size() > 3)
        inputRead.setMaxTime( atoi( fline[3].c_str() ));
    }
//    else if ( fline[1] == "mba" )
//    {
////      inputRead.setNBodies( atoi( fline[2].c_str() ) ) ;
//      inputRead.setDistance( atof( fline[2].c_str() ) ) ;
//    }
  } else if (keyword == "omp")
  {
    cout << "OMP command found" << endl;
    inputRead.setOMP(atoi(fline[1].c_str()));
  } else if (keyword == "pbc")
  {
    cout << "PBC command found" << endl;
    inputRead.setPBCT( atoi(fline[1].c_str()) );
    if ( inputRead.getPBCs() > 0 )
      inputRead.setBoxl(atof(fline[2].c_str()));
//  } else if (keyword == "cutoffs" ) {
//      inputRead.setCutoffInter(atof(fline[1].c_str()));
//      inputRead.setCutoffInteract(atof(fline[2].c_str()));
//      inputRead.setCutoffIntra(atof(fline[3].c_str()));
  } else if (keyword == "salt")
  {
    cout << "Salt command found" << endl;
    inputRead.setSaltCon( atof(fline[1].c_str()) );
  } else if (keyword == "temp")
  {
    cout << "Temperature command found" << endl;
    inputRead.setTemp( atof(fline[1].c_str()) );
  } else if (keyword == "idiel")
  {
    cout << "Interior dielectric command found" << endl;
    inputRead.setIDiel( atof(fline[1].c_str()) );
  } else if (keyword == "sdiel")
  {
    cout << "Solvent dielectric command found" << endl;
    inputRead.setSDiel( atof(fline[1].c_str()) );
  } else if (keyword == "attypes")
  {
    cout << "Atom Types command found" << endl;
    inputRead.setNType( atoi(fline[1].c_str()) );
    inputRead.resizeVecs();
    cout << "done with attypes" << endl;
  } else if (keyword == "type")
  {
    cout << "Type def command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > inputRead.getType()-1)
      return;
    if (fline.size() > 2)
      inputRead.setTypeNCount( typeNo, atoi(fline[2].c_str()) );
    if (fline.size() > 3)
    {
      inputRead.setTypeNDef( typeNo, fline[3].c_str() );
      if (inputRead.getTypeNDef(typeNo) == "move")
        inputRead.setTypeNDtr( typeNo, atof(fline[4].c_str()));
    }
  } else if (keyword == "pqr")
  {
    cout << "PQR command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > inputRead.getType()-1)
      return;
    inputRead.setTypeNPQR( typeNo, fline[2].c_str() );
//  } else if (keyword == "imat")
//  {
//    cout << "IMAT command found" << endl;
//    int typeNo = atoi(fline[1].c_str())-1;
//    if (typeNo > inputRead.getType()-1)
//      return;
//    inputRead.setTypeNImat( typeNo, fline[2] );
//  } else if (keyword == "spolDir")
//  {
//    cout << "SelfpolDir command found" << endl;
//    int typeNo = atoi(fline[1].c_str())-1;
//    if (typeNo > inputRead.getType()-1)
//      return;
//    inputRead.setTypeNSpolDir( typeNo, fline[2] );
//  } else if (keyword == "spolExt")
//  {
//    cout << "SpolName command found" << endl;
//    int typeNo = atoi(fline[1].c_str())-1;
//    if (typeNo > inputRead.getType()-1)
//      return;
//    inputRead.setTypeNExp( typeNo, fline[2] );
  } else if (keyword == "xyz")
  {
    cout << "XYZ command found" << endl;
    int typeNo = atoi(fline[1].c_str())-1;
    if (typeNo > inputRead.getType()-1)
      return;
    inputRead.setTypeNXYZ( typeNo, fline[2] );
//  } else if (keyword == "f_ext")
//  {
//    cout << "External force command found" << endl;
//    inputRead.setFExt( atof(fline[1].c_str()), atof(fline[2].c_str()));
//  } else if (keyword == "scale")
//  {
//    cout << "Scale command found" << endl;
//    inputRead.setScale( atof(fline[1].c_str()));
//  } else if (keyword == "distance")
//  {
//    cout << "Distance command found" << endl;
//    inputRead.setMolDist( atof(fline[1].c_str()));
//  } else if (keyword == "layer")
//  {
//    cout << "Layer command found" << endl;
//    inputRead.setLayer( atoi(fline[1].c_str()));
  } else if (keyword == "random")
  {
    cout << "RNG Seed command found" << endl;
    inputRead.setRand( atoi(fline[1].c_str()) );
  } else
    cout << "Keyword not found, read in as " << fline[0] << endl;
} // end findKeyword

/*!#########################################################*/
/*#########################################################*/
// readInputFile
// Read an input file and return parameters for various words
/*#########################################################*/
/*#########################################################*/
void readInputFile(char * fname, CSetup & readIn)
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
    findLines(inputLine, readIn);
    getline(fin, inputLine);
  }
} // end readInputFile

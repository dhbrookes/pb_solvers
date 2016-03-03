//
//  InputFile.cpp
//  makespheres
//
//  Created by David Brookes on 8/25/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#include "InputFile.h"

double DEFAULT_TOLSP = 5;
double DEFAULT_TOLIP = 1;

/*#########################################################*/
//  Constructor for BaseInputFile
/*#########################################################*/
BaseInputFile::BaseInputFile(string path)
:path_(path)
{
}

/*#########################################################*/
//  read
// function to read a file, search for keywords
/*#########################################################*/
void BaseInputFile::read()
{
  
  FILE* file = fopen(path_.c_str(), "r"); /* should check the result */
  cout << "\nReading input file: " << path_.c_str() << endl;
  
  char line[256];
  char delim = ' ';
  vector<string> tokens;
  string key, val, s;
  char* line_str = (char*)malloc(256);

  while (fgets(line, sizeof(line), file))
  {
    
    line_str = &strtok(line, "\n")[0]; // remove new line
    if (line_str == NULL) continue;
    vector<string> tokens = line_tokenize(line_str, delim);
    key = tokens[0];
    val = tokens[1];
    for(int i = 0; i < KEYWORDS.size(); i++)
    {
      s = KEYWORDS[i];
      if (key == s) search_keywords(s, val);
    }
  }
  fclose(file);
} // end BaseInputFile::read

/*#########################################################*/
//  search_keywords
// Reads a string and searches for keyword in it
/*#########################################################*/
void AllInputFile::search_keywords(string s, string val)
{
    if (s == KEYWORDS[1])
    {
      cout << "PDB command found" << endl;
      pdb_ = strip_quotes(val);
      in_sat = true;
    }
    else if (s == KEYWORDS[2])
    {
      cout << "PQR command found" << endl;
      pqr_ = strip_quotes(val);
      in_sat = true;
      makePQR_ = false;
    }
    else if (s == KEYWORDS[3])
    {
      cout << "pH command found" << endl;
      pH_ = strtod(val.c_str(), NULL);
    }
    else if (s == KEYWORDS[4])
    {
      cout << "outfiledir command found" << endl;
      outfiledir_ = strip_quotes(val);
      home_sat = true;
    }
    else if (s == KEYWORDS[5])
    {
      cout << "tempfiledir command found" << endl;
      tempfiledir_ = strip_quotes(val);
      home_sat = true;
    }
} // end AllInputFile::search_keywords

/*#########################################################*/
// AllInputFile constructor
// identifies list of accepted keywords, and their default
// values
/*#########################################################*/
AllInputFile::AllInputFile(string path)
:BaseInputFile(path)
{
  KEYWORDS.push_back("homedir");//"HOMEDIR");
  KEYWORDS.push_back("pdb");//"PDB");
  KEYWORDS.push_back("pqr");//"PQR");
  KEYWORDS.push_back("ph");
  KEYWORDS.push_back("outfiledir");//"OUTFILEDIR");
  KEYWORDS.push_back("tempfiledir");
  homedir_ = "../../src";
  pH_=8; //default
  home_sat = true;
  in_sat = false;
  makePQR_ = true;
  tempfiledir_ = "None";
  outfiledir_ = "None";
  
  read();
  
  if (tempfiledir_ == "None") tempfiledir_ = homedir_;
  if (outfiledir_ == "None") outfiledir_ = homedir_;
  else if (!in_sat) 
    throw NotEnoughInputsException("need an input .pdb or .pqr file with the \
    keywords 'pdb' or 'pqr'");
  
}// end AllInputFile::AllInputFile


//
//  InputFile.h
//  makespheres
//
//  Created by David Brookes on 8/25/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#ifndef __makespheres__InputFile__
#define __makespheres__InputFile__

#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cstring>
#include "readutil.h"

using namespace std;

class BaseInputFile
{
protected:
  class NotEnoughInputsException : public std::exception
  {
    public:
    string message_;
    NotEnoughInputsException(string message=" ")
    :message_(message)  {  }
  
    virtual const char* what() const throw()
    {
      string full = "Not enough inputs: " + message_;
      return full.c_str();
    }
    virtual ~NotEnoughInputsException() throw() {}
  };
  
  string path_;
  vector<string> KEYWORDS;
  
  void read();
  virtual void search_keywords(string s, string val)=0;
  
  BaseInputFile(string path);
  
  const std::string& get_path() const { return path_; }
}; // end BaseInputFile


class AllInputFile : public BaseInputFile
{
    
protected:
  
  string homedir_; //path to $PBAMHOME
  string pdb_;
  string pqr_; //if the pqr is pre-made
  bool makePQR_; //true if pdb2pqr needs to be run
  double pH_;
  string tempfiledir_; //dir to store intermed files from progs. Def: homedir_
  string outfiledir_;
  
  bool home_sat; //homedir condition satisfied
  bool in_sat; //input file condition satisfied
  
public:
  
  AllInputFile(string path);
  void search_keywords(string s, string val);
  void printVals();
  
  const string& get_homedir() const           { return homedir_ ;         }
  const string& get_pdb() const               { return pdb_;              }
  const string& get_pqr() const               { return pqr_;              }
  const bool get_make_pqr() const             { return makePQR_;          }
  const double get_ph() const                 { return pH_;               }
  const string& get_tempfile_dir() const      { return tempfiledir_;      }
  
}; // end AllInputFile

#endif /* defined(__makespheres__InputFile__) */


//
//  pdb2pqr_wrap.cpp
//  all
//
//  Created by David Brookes on 8/24/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#include "pdb2pqr_wrap.h"
using namespace std;


PDB2PQR_WRAP::PDB2PQR_WRAP(const string& script_path,
                           const string& in_path,
                           const string& out_path,
                           const double pH,
                           const FF_OPTION ff
                           )
:SCRIPT_PATH(script_path), pH_(pH), inPath_(in_path), 
ff_(ff), outPath_(out_path)
{
  string ff_str;
  
  if (ff_ == AMBER)
  {
      ff_str="AMBER";
  }
  else
  {
      ff_str="AMBER";
  }
  
  ostringstream strs;
  strs << pH_;
  std::string str1 = strs.str();
  strs.clear();
  
  command_ = "python " + SCRIPT_PATH + " --chain";
  command_ += " --with-ph=" + str1;
  command_ += " --ff=" + ff_str;
  command_ += " " + inPath_;
  command_ += " " + outPath_;
}


const string PDB2PQR_WRAP::run() const
{
  system(command_.c_str());
  return outPath_;
}

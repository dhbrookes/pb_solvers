//
//  pdb2pqr_wrap.h
//  all
//
//  Created by David Brookes on 8/24/15.
//  Copyright (c) 2015 David Brookes. All rights reserved.
//

#ifndef __all__pdb2pqr_wrap__
#define __all__pdb2pqr_wrap__

#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>

using namespace std;

enum FF_OPTION { AMBER };

class PDB2PQR_WRAP
{
protected:
    
    string SCRIPT_PATH;
    double pH_;
    FF_OPTION ff_;
    string inPath_;
    string outPath_;
    string command_;

public:
    PDB2PQR_WRAP(const string& script_path,
                 const string& in_path,
                 const string& out_path,
                 const double pH=8,
                 const FF_OPTION ff=AMBER);
    
    const string run() const;
};


#endif /* defined(__all__pdb2pqr_wrap__) */

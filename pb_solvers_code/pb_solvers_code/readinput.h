#ifndef READINP_H
#define READINP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include "readutil.h"
#include "setup.h"
#include "util.h"


vector<string> split(string str, char delimiter);
void findLines(string fline, CSetup & inputRead);
void findKeyword(vector<string> fline, CSetup & inputRead);
void readInputFile(char * fname, CSetup & readIn);

#endif


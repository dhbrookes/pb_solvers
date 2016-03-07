#ifndef READINP_H
#define READINP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include "setup.h"


vector<string> split(string str, char delimiter);
void findLines(string fline, CSetup & inputRead);
void findKeyword(vector<string> fline, CSetup & inputRead);
void readInputFile(string fname, CSetup & readIn);

#endif


//
//  readutil.h
//  pb_solvers_code
//
//  Created by David Brookes on 11/18/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef readutil_h
#define readutil_h

#include "util.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

#include <string>

using namespace std;

class CouldNotReadException: public exception
{
protected:
  string path_;
  
public:
  CouldNotReadException(string path)
  :path_(path)
  {
  }
  
  virtual const char* what() const throw()
  {
    string ss;
    ss = "Could not read: " + path_;
    return ss.c_str();
  }
};

class NotEnoughCoordsException: public exception
{
protected:
  int infile_;
  int needed_;
  
public:
  NotEnoughCoordsException( int inFile, int needed)
  :infile_(inFile), needed_(needed)
  {
  }
  
  virtual const char* what() const throw()
  {
    string ss;
    ss = "File has "+to_string(infile_)+" lines, need "+to_string(needed_);
    return ss.c_str();
  }
  
};

/*
 Class for reading and storing the info in .pqr files
 */
class PQRFile
{
protected:
  string path_;
  int M_;  // number of charges
  vector<double> charges_;
  vector<double> atomRadii_;
  vector<Pt> atomCenters_;
  
  bool cg_;
  vector<Pt> cgCenters_;  // coarse grain centers
  vector<double> cgRadii_;
  
  void read()
  {
    ifstream fin(path_.c_str());
    if (!fin.is_open()) throw CouldNotReadException(path_);
    
    char buf[600];
    fin.getline(buf, 599);
    Pt cent( 0.0, 0.0, 0.0);
    int iCen = 17;
    int iCoord = 31;
    while (!fin.eof())
    {
      double x, y, z, c, r;
      if (strncmp(&(buf[0]),"ATOM",4) == 0)
      {
        sscanf(&(buf[iCoord]), "%lf %lf %lf %lf %lf", &x, &y, &z, &c, &r);
        // read in as centers that specifies dielectric boundary
        if (strncmp(&(buf[iCen]),"CEN",3) == 0)
        {
          cg_ = true;
          cgRadii_.push_back(r);
          cgCenters_.push_back(Pt(x,y,z));
        }
        
        // read in as atoms
        else
        {
          atomCenters_.push_back(Pt(x,y,z));
          charges_.push_back(c);
          atomRadii_.push_back(r);
          cent = cent + Pt(x,y,z);
        }
      }
      fin.getline(buf,99);
    }
    M_ = int(atomCenters_.size());
    if (cg_ == false) cgCenters_.push_back( cent * (1./(double) M_));
  }
  
public:
  
  PQRFile(string path, int approx_size=1000)
  :path_(path), cg_(false)
  {
    charges_.reserve(approx_size);
    atomRadii_.reserve(approx_size);
    atomCenters_.reserve(approx_size);
    cgCenters_.reserve(10);
    cgRadii_.reserve(10);
    read();
  }
  
  const string get_path() const             { return path_; }
  vector<Pt> get_atom_pts() const           { return atomCenters_; }
  const vector<double> get_charges() const  { return charges_; }
  const vector<double> get_radii() const    { return atomRadii_; }
  const int get_M() const                   { return M_; }
  const vector<double> get_cg_radii() const { return cgRadii_; }
  const vector<Pt> get_cg_centers() const   { return cgCenters_; }
  const bool get_cg() const                 { return cg_; }
  
};


/*
 Class for reading and storing the info in .xyz files
 */
class XYZFile
{
protected:
  int nmols_;  // number of molecules
  vector<Pt> pts_;
  string path_;
  
public:
  XYZFile(string path, int nmols)
  :path_(path), nmols_(nmols)
  {
    pts_.reserve(nmols);
    read();
  }
  
  void read()
  {
    int ctr = 0;
    string s;
    double x,y,z;
    ifstream fin (path_.c_str());
    if (!fin.is_open()) throw CouldNotReadException(path_);
    
    while ( getline( fin, s ) && (ctr < nmols_))
    {
      istringstream ss( s );
      ss >> x >> y >> z;
      pts_.push_back(Pt(x, y, z));
      ctr++;
    }
    if (ctr != nmols_)  throw NotEnoughCoordsException(ctr, nmols_);
  }
  
  const string get_path() const     { return path_; }
  const int get_nmols() const       { return nmols_; }
  vector<Pt> get_pts() const        { return pts_; }
};
    

#endif /* readutil_h */

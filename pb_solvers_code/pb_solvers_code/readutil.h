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
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

//vector<string> line_tokenize(string line, char delim)
//{
//  vector<string> v;
//  istringstream buf(line);
//  for(string token; getline(buf, token, delim); )
//  {
//    v.push_back(token);
//  }
//  return v;
//}


//string strip_quotes(string s)
//{
//  s.erase(remove( s.begin(), s.end(), '\"' ), s.end());
//  return s;
//}

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
    ostringstream ss;
    ss << "Could not read: " << path_ << endl;
    return ss.str().c_str();
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
  vector<Pt> pts_;
  vector<double> charges_;
  vector<double> atomRadii_;
  vector<Pt> atomCenters_;
  
  bool cg_;
  vector<Pt> cgCenters_;  // coarse grain centers
  vector<double> cgRadii_;
  
  void read()
  {
    ifstream fin(path_.c_str());
    if (!fin.is_open())
    {
      throw CouldNotReadException(path_);
    }
    char buf[100];
    fin.getline(buf, 99);
    int iCen = 18;
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
        }
      }
      
      fin.getline(buf,99);
    }
    M_ = int(pts_.size());
  }
  
public:
  
  PQRFile(string path, int aprox_size=100)
  :path_(path), pts_(aprox_size), charges_(aprox_size), atomRadii_(aprox_size),
    atomCenters_(aprox_size), cgCenters_(10), cgRadii_(10), cg_(false)
  {
    read();
  }
  
  const string get_path() const             { return path_; }
  const vector<Pt> get_pts() const          { return pts_; }
  const vector<double> get_charges() const  { return charges_; }
  const vector<double> get_radii() const    { return atomRadii_; }
  const vector<Pt> get_centers() const      { return atomCenters_; }
  const int get_M() const                   { return M_; }
  const vector<double> get_cg_radii() const { return cgRadii_; }
  const vector<Pt> get_cg_centers() const    { return cgCenters_; }
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
    string s;
    double x,y,z;
    ifstream fin (path_.c_str());
    if (!fin.is_open())
    {
      throw CouldNotReadException(path_);
    }
    
    for (int i = 0; i < nmols_; i++)
    {
      getline( fin, s );
      istringstream ss( s );
      ss >> x >> y >> z;
      pts_.push_back(Pt(x, y, z));
    }
  }
  
  const string get_path() const     { return path_; }
  const int get_nmols() const       { return nmols_; }
  const vector<Pt> get_pts() const  { return pts_; }
};
    

#endif /* readutil_h */

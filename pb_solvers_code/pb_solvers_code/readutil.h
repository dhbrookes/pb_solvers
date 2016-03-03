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

vector<string> line_tokenize(string line, char delim)
{
  vector<string> v;
  istringstream buf(line);
  for(string token; getline(buf, token, delim); )
  {
    v.push_back(token);
  }
  return v;
}


string strip_quotes(string s)
{
  s.erase(remove( s.begin(), s.end(), '\"' ), s.end());
  return s;
}

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
  vector<Pt> pts_;
  vector<double> charges_;
  vector<double> radii_;
  vector<Pt> centers_;
  
public:
  PQRFile(string path, int aprox_size=100)
  :path_(path)
  {
    pts_.reserve(aprox_size);
    charges_.reserve(aprox_size);
    radii_.reserve(aprox_size);
    centers_.reserve(aprox_size);
    
    read();
  }
  
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
      double x,y,z,c,r;
      sscanf(&(buf[iCoord]), "%lf %lf %lf %lf %lf", &x, &y, &z, &c, &r);
      // read in as centers that specifies dielectric boundary:
      if (strncmp(&(buf[iCen]),"CEN",3) == 0)
      {
        radii_.push_back(r);
        centers_.push_back(Pt(x,y,z));
        
      }
      
      // read in as atoms
      else
      {
        pts_.push_back(Pt(x,y,z));
        charges_.push_back(c);
      }
      fin.getline(buf,99);
    }
  }
  
  const string get_path() const             { return path_; }
  const vector<Pt> get_pts() const          { return pts_; }
  const vector<double> get_charges() const  { return charges_; }
  const vector<double> get_radii() const    { return radii_; }
  const vector<Pt> get_centers() const      { return centers_; }
  
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

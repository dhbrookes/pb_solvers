//
//  Solver.h
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef Solver_h
#define Solver_h

#include <stdio.h>
#include "util.h"
#include <memory>
#include "System.h"
#include "SHCalc.h"

/*
 Pre-computed values representing fixed charges within each sphere. Each object
 of this class refers to one molecule. This is then a vector of matrices where
 the vector index, k, loops through each coarse-grained sphere in the molecule.
 The indices of the inner matrices are then n and m, which loop over the number
 of poles in the system. See eq. 8a in Yap 2010 for more info
 */
class EMatrix
{
protected:
  vector<MyMatrix<cmplx> > E_;
  int p_;
  
  void set_E_knm(int k, int n, int m, cmplx val) { E_[k].set_val(n, m+p_, val); }
  
public:
  EMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc, int p, double eps_in);
  
  cmplx get_E_knm(int k, int n, int m) { return E_[k](n, m+p_); }
  MyMatrix<cmplx> get_E_k(int k)     { return E_[k]; }
};

/*
 Pre-computed values representing fixed charges outside each sphere. Each object
 of this class refers to one molecule. This is then a vector of matrices where
 the vector index, k, loops through each coarse-grained sphere in the molecule.
 The indices of the inner matrices are then n and m, which loop over the number
 of poles in the system. See eq. 8b in Yap 2010 for more info
 */
class LEMatrix
{
protected:
  vector<MyMatrix<cmplx> > LE_;
  int p_;
  
  void set_LE_knm(int k, int n,
                  int m, cmplx val) { LE_[k].set_val(n, m+p_, val); }
  
public:
  LEMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc, int p, double eps_in);
  
  cmplx get_LE_knm(int k, int n, int m) { return LE_[k](n, m+p_); }
  MyMatrix<cmplx> get_LE_k(int k)     { return LE_[k]; }
};


/*
 Class for pre-computing values of surface integral matrices I_E. Each object 
 of this class refers to one molecule. See equation 21 in Yap 2010 for more
 info.
 */
class IEMatrix
{
protected:
  
  // indices in order are k, (n, m), (l, s)
  vector<MatOfMats<cmplx>::type > IE_;
  int p_;
  
  void set_IE_k_nm_ls(int k, int n, int m, int l, int s, cmplx val)
  {
    IE_[k](n, m+p_).set_val(l, s+p_, val);
  }
  
  cmplx get_IE_k_nm_ls(int k, int n, int m, int l, int s)
  {
    return IE_[k](n, m+p_)(l, s+p_);
  }
  
  
  /*
   Make a uniform grid of points on the surface of the sphere. Based on the 
   Fibonacci sphere algorithm
   */
  vector<Pt> make_uniform_sph_grid(int num, double r);
  
public:
  
  IEMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc, int p);
  
};

#endif /* Solver_h */

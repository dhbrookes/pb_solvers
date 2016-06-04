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
#include "BesselCalc.h"
#include "ReExpCalc.h"

/*
 Base class for matrices that will be re-expanded
 */
class ComplexMoleculeMatrix
{
protected:
  vector<MyMatrix<cmplx> > mat_;
  int p_;
  int I_;
  
  void set_mat_knm(int k, int n, int m, cmplx val) { mat_[k].set_val(n, m+p_, val); }
    
public:
  ComplexMoleculeMatrix(int I, int ns, int p);
  
  cmplx get_mat_knm(int k, int n, int m) { return mat_[k](n, m+p_); }
  const int get_I() const   { return I_; }
  const int get_p() const   { return p_; }
    
};

/*
 Pre-computed values representing fixed charges within each sphere. Each object
 of this class refers to one molecule. This is then a vector of matrices where
 the vector index, k, loops through each coarse-grained sphere in the molecule.
 The indices of the inner matrices are then n and m, which loop over the number
 of poles in the system. See eq. 8a in Yap 2010 for more info
 */
class EMatrix: public ComplexMoleculeMatrix
{
public:
  EMatrix(int I, int ns, int p);
  virtual void calc_vals(Molecule mol, shared_ptr<SHCalc> sh_calc,
                         double eps_in);
  
};

/*
 Pre-computed values representing fixed charges outside each sphere. Each object
 of this class refers to one molecule. This is then a vector of matrices where
 the vector index, k, loops through each coarse-grained sphere in the molecule.
 The indices of the inner matrices are then n and m, which loop over the number
 of poles in the system. See eq. 8b in Yap 2010 for more info
 */
class LEMatrix : public ComplexMoleculeMatrix
{
public:
  LEMatrix(int I, int ns, int p);
  void calc_vals(Molecule mol, shared_ptr<SHCalc> sh_calc,
                 double eps_in);
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
  int I_;
  
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
  IEMatrix(int I, int ns, int p);

  IEMatrix(int I, Molecule mol, shared_ptr<SHCalc> sh_calc, int p);
  
  void calc_vals(Molecule mol, shared_ptr<SHCalc> sh_calc);
};


/*
 Re-expansion coefficients
 */
class TMatrix
{
protected:
  int p_;
  double kappa_;
  vector<shared_ptr<ReExpCoeffs> > T_;
  // maps (I,k), (J,l) indices to a single index in T. If this returns -1
  // then the spheres are overlapping or less than 5A away and numerical
  // re-expansion is required
  map<vector<int>, int> idxMap_;
  
public:
  
  TMatrix(int p, shared_ptr<System> _sys, shared_ptr<SHCalc> _shcalc,
          shared_ptr<Constants> _consts, shared_ptr<BesselCalc> _besselcalc,
          shared_ptr<ReExpCoeffsConstants> _reexpconsts);
  
  // if these spheres can be re-expanded analytically, return true
  bool is_analytic(int I, int k, int J, int l)
  {
    if (idxMap_[{I, k, J, l}] == -1 ) return false;
    else return true;
  }
  
  // get re-expansion of sphere (I, k) with respect to (J, l)
  shared_ptr<ReExpCoeffs> get_T_Ik_Jl(int I, int k, int J, int l)
  {
    return T_[idxMap_[{I, k, J, l}]];
  }
  
  void update_vals(shared_ptr<System> _sys, shared_ptr<SHCalc> _shcalc,
                   shared_ptr<BesselCalc> _besselcalc,
                   shared_ptr<ReExpCoeffsConstants> _reexpconsts);
  
  
  /*
   Re-expand a matrix X with respect to T(I,k)(J,l)
  */
  MyMatrix<cmplx> re_expand(int I, int k, int J, int l, MyMatrix<cmplx> X);
  
};

class HMatrix;
class FMatrix;

/*
 Equation 8c
 */
class LFMatrix : public ComplexMoleculeMatrix
{
public:
  LFMatrix(int I, int ns, int p);
  
};

/*
 Equation 10b
 */
class LHMatrix : public ComplexMoleculeMatrix
{
public:
  LHMatrix(int I, int ns, int p);
};

/*
 Equation 10c
 */
class LHNMatrix : public ComplexMoleculeMatrix
{
public:
  LHNMatrix(int I, int ns, int p);
  
};

/*
 Equation 14a
 */
class XHMatrix : public ComplexMoleculeMatrix
{
public:
  XHMatrix(int I, int ns, int p);
  
  void calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                 shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                 shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                 shared_ptr<LHNMatrix> LHN);

};

/*
 Equation 14b
 */
class XFMatrix : public ComplexMoleculeMatrix
{
public:
  XFMatrix(int I, int ns, int p);
  
  void calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                 shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                 shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                 shared_ptr<LHNMatrix> LHN, double eps_in,
                 double eps_out, double kappa);
};


#endif /* Solver_h */




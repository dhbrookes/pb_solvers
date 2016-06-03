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
  int I_;
  
  void set_E_knm(int k, int n, int m, cmplx val) { E_[k].set_val(n, m+p_, val); }
  
public:
  //i is index of molecule, ns is number of coarse grained spheres in molecule
  EMatrix(int I, int ns, int p);
  EMatrix(int I, Molecule mol, shared_ptr<SHCalc> sh_calc, int p, double eps_in);
  
  void calc_vals(Molecule mol, shared_ptr<SHCalc> sh_calc, double eps_in);
  
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
  int I_;
  
  void set_LE_knm(int k, int n,
                  int m, cmplx val) { LE_[k].set_val(n, m+p_, val); }
  
public:
  LEMatrix(int I, int ns, int p);
  LEMatrix(int I, Molecule mol, shared_ptr<SHCalc> sh_calc, int p, double eps_in);
  
  void calc_vals(Molecule mol, shared_ptr<SHCalc> sh_calc, double eps_in);
  
  cmplx get_LE_knm(int k, int n, int m) { return LE_[k](n, m+p_); }
  MyMatrix<cmplx> get_LE_k(int k)       { return LE_[k]; }
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
  bool is_analytic(int I, int k, int J, int l);
  
  // get re-expansion of sphere (I, k) with respect to (J, l)
  shared_ptr<ReExpCoeffs> get_T_Ik_Jl(int I, int k, int J, int l);
  
  void update_vals(shared_ptr<System> _sys, shared_ptr<SHCalc> _shcalc,
                   shared_ptr<BesselCalc> _besselcalc,
                   shared_ptr<ReExpCoeffsConstants> _reexpconsts);
  
};


class HMatrix;
class FMatrix;

/*
 Equation 8c
 */
class LFMatrix
{
protected:
  vector<MyMatrix<cmplx> > LF_;
  int p_;
  int I_;
  
  void set_LF_knm(int k, int n,
                  int m, cmplx val) { LF_[k].set_val(n, m+p_, val); }
  
public:
  
  LFMatrix(int I, int ns, int p);
  
  cmplx get_LF_knm(int k, int n, int m)   { return LF_[k](n, m+p_); }
  MyMatrix<cmplx> get_LF_k(int k)         { return LF_[k]; }
  
};

/*
 Equation 10b
 */
class LHMatrix
{
protected:
  vector<MyMatrix<cmplx> > LH_;
  int p_;
  int I_;
  
  void set_LH_knm(int k, int n,
                  int m, cmplx val) { LH_[k].set_val(n, m+p_, val); }
  
public:
  
  LHMatrix(int I, int ns, int p);
  
  cmplx get_LH_knm(int k, int n, int m)   { return LH_[k](n, m+p_); }
  MyMatrix<cmplx> get_LH_k(int k)         { return LH_[k]; }
  
};

/*
 Equation 10c
 */
class LHNMatrix
{
protected:
  vector<MyMatrix<cmplx> > LHN_;
  int p_;
  int I_;
  
  void set_LHN_knm(int k, int n,
                  int m, cmplx val) { LHN_[k].set_val(n, m+p_, val); }
  
public:
  
  LHNMatrix(int I, int ns, int p);
  
  cmplx get_LHN_knm(int k, int n, int m)   { return LHN_[k](n, m+p_); }
  MyMatrix<cmplx> get_LHN_k(int k)         { return LHN_[k]; }
  
};

/*
 Equation 14a
 */
class XHMatrix
{
protected:
  vector<MyMatrix<cmplx> > XH_;
  int p_;
  int I_;
  
  void set_XH_knm(int k, int n,
                   int m, cmplx val) { XH_[k].set_val(n, m+p_, val); }
  
public:
  XHMatrix(int I, int ns, int p);
  
  void calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                 shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                 shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                 shared_ptr<LHNMatrix> LHN);
  
  cmplx get_XH_knm(int k, int n, int m)   { return XH_[k](n, m+p_); }
};


/*
 Equation 14b
 */
class XFMatrix
{
protected:
  vector<MyMatrix<cmplx> > XF_;
  int p_;
  int I_;
  
  void set_XF_knm(int k, int n,
                  int m, cmplx val) { XF_[k].set_val(n, m+p_, val); }
  
public:
  XFMatrix(int I, int ns, int p);
  
  void calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                 shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                 shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                 shared_ptr<LHNMatrix> LHN, double eps_in,
                 double eps_out, double kappa);
  
  cmplx get_XF_knm(int k, int n, int m)   { return XF_[k](n, m+p_); }
};


#endif /* Solver_h */




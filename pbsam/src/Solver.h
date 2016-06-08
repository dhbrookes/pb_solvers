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
#include <memory>
#include "ReExpCalc.h"

// use Rakhmanov method
vector<Pt> make_uniform_sph_grid(int m_grid, double r);

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
  MyMatrix<cmplx> get_mat_k(int k)       { return mat_[k]; }
  const int get_I() const   { return I_; }
  const int get_p() const   { return p_; }
  
  void reset_mat();
    
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

  /*
   Make a uniform grid of points on the surface of the sphere.
   */
//  vector<Pt> make_uniform_sph_grid(int num, double r);
  
public:
  IEMatrix(int I, int ns, int p);

  IEMatrix(int I, Molecule mol, shared_ptr<SHCalc> sh_calc, int p);
  
  cmplx get_IE_k_nm_ls(int k, int n, int m, int l, int s)
  {
    return IE_[k](n, m+p_)(l, s+p_);
  }
  
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
  shared_ptr<SHCalc>  _shCalc_;
  shared_ptr<BesselCalc>  _besselCalc_;
  
  int Nmol_;
  vector<int> Nsi_; // number of spheres in each molecule
  
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
  
  int get_nmol() const { return Nmol_; }
  
  int get_nsi(int i)   { return Nsi_[i]; }
  
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
  
  void calc_vals(shared_ptr<TMatrix> T, shared_ptr<FMatrix> F,
                 shared_ptr<SHCalc> shcalc, shared_ptr<System> sys);
  
  // analytic re-expansion (Equation 27a)
  MyMatrix<cmplx> analytic_reex(int I, int k, int j,
                                shared_ptr<FMatrix> F,
                                shared_ptr<SHCalc> shcalc,
                                shared_ptr<System> sys, int Mp=-1);
  
  /*
   Equation 15a. For analytic re expansion
   */
  cmplx make_fb_Ij(int I, int j, Pt rb,
                   shared_ptr<FMatrix> F,
                   shared_ptr<SHCalc> shcalc);
};

/*
 Equation 10b
 */
class LHMatrix : public ComplexMoleculeMatrix
{
protected:
  double kappa_;
  
public:
  LHMatrix(int I, int ns, int p, double kappa);
  
  void calc_vals(shared_ptr<TMatrix> T, shared_ptr<HMatrix> H,
                 shared_ptr<SHCalc> shcalc, shared_ptr<System> sys,
                 shared_ptr<BesselCalc> bcalc, int Mp=-1);
  
  // analytic re-expansion (Equation 27b)
  MyMatrix<cmplx> analytic_reex(int I, int k, int j,
                                shared_ptr<HMatrix> H,
                                shared_ptr<SHCalc> shcalc,
                                shared_ptr<System> sys,
                                shared_ptr<BesselCalc> bcalc,
                                int Mp=-1);
  
  /*
   Equation 15b. For analytic re expansion
   */
  cmplx make_hb_Ij(int I, int j, Pt rb,
                   shared_ptr<HMatrix> H,
                   shared_ptr<SHCalc> shcalc,
                   shared_ptr<BesselCalc> bcalc);
  
};

/*
 Equation 10c
 */
class LHNMatrix : public ComplexMoleculeMatrix
{
public:
  LHNMatrix(int I, int ns, int p);
  
  void calc_vals(shared_ptr<TMatrix> T, vector<shared_ptr<HMatrix> > H);
  
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
protected:
  double eps_;
  
public:
  XFMatrix(int I, int ns, int p, double eps_in, double eps_out);
  
  void calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                 shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                 shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                 shared_ptr<LHNMatrix> LHN, double kappa);
  
  const double get_eps() const { return eps_; }
};


class HMatrix: public ComplexMoleculeMatrix
{
protected:
  double kappa_;
  
public:
  HMatrix(int I, int ns, int p, double kappa);
  
  void calc_vals(Molecule mol, shared_ptr<HMatrix> prev,
                 shared_ptr<XHMatrix> XH,
                 shared_ptr<FMatrix> F,
                 shared_ptr<IEMatrix> IE,
                 shared_ptr<BesselCalc> bcalc);
  
};


class FMatrix: public ComplexMoleculeMatrix
{
protected:
  double kappa_;
  
public:
  FMatrix(int I, int ns, int p, double kappa);
  
  void calc_vals(Molecule mol, shared_ptr<FMatrix> prev,
                 shared_ptr<XFMatrix> XF,
                 shared_ptr<HMatrix> H,
                 shared_ptr<IEMatrix> IE,
                 shared_ptr<BesselCalc> bcalc);
  
};

#endif /* Solver_h */




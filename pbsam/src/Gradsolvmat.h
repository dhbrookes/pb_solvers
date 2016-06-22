//
//  Gradsolvmat.h
//  pbsam_xcode
//
//  Created by David Brookes on 6/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef Gradsolvmat_h
#define Gradsolvmat_h

#include <stdio.h>
#include "MyMatrix.h"
#include "util.h"
#include "Solvmat.h"

/*
 References:
 [1] Yap, E., Head-Gordon, T. 2010. JCTC
 [2] Yap, E., Head-Gordon, T. 2013. JCTC
 */

/*
 Base class for gradients of ComplexMoleculeMatrix objects
 */
class GradCmplxMolMat
{
protected:
  int wrt_;  // with respect to
  vector<MyMatrix<Ptx> > mat_;  // each gradient has three components (use Pt)
  int p_;
  int I_;
  
  
public:
  GradCmplxMolMat(int I, int wrt, int ns, int p)
  :wrt_(wrt), p_(p), I_(I), mat_(ns, MyMatrix<Ptx> (p, 2*p+1))
  {
  }
  
  const int get_wrt() const   { return wrt_; }
  const int get_I() const   { return I_; }
  const int get_p() const   { return p_; }
  const int get_ns() const  { return (int) mat_.size(); }
  Ptx get_mat_knm(int k, int n, int m) { return mat_[k](n, m+p_); }
  void set_mat_knm(int k, int n, int m, Ptx val)
  { mat_[k].set_val(n, m+p_, val); }
  
  MyMatrix<Ptx> get_mat_k(int k) const { return mat_[k]; }
  
  void reset_mat();
};

// gradient matrices:
class GradFMatrix;
class GradHMatrix;
class GradWHMatrix;
class GradWFMatrix;
class GradLHNMatrix;
class GradLHMatrix;
class GradLFMatrix;

/*
 Eq. 12a [2]
 */
class GradWFMatrix : public GradCmplxMolMat
{
protected:
  double eps_;
  double kappa_;
  
  
public:
  GradWFMatrix(int I, int wrt, int ns, int p,
               double eps_in, double eps_out,
               double kappa);
  
  // calculate values for all spheres in this molecule
  void calc_all_vals(shared_ptr<Molecule> mol,
                     shared_ptr<PrecalcBessel> bcalc,
                     shared_ptr<GradHMatrix> dH,
                     shared_ptr<GradFMatrix> dF,
                     shared_ptr<GradLHMatrix> dLH,
                     shared_ptr<GradLHNMatrix> dLHN,
                     shared_ptr<GradLFMatrix> dLF);
  
  // calcualte values for one sphere in this molecule
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  shared_ptr<PrecalcBessel> bcalc,
                  shared_ptr<GradHMatrix> dH,
                  shared_ptr<GradFMatrix> dF,
                  shared_ptr<GradLHMatrix> dLH,
                  shared_ptr<GradLHNMatrix> dLHN,
                  shared_ptr<GradLFMatrix> dLF);
  
  const double get_eps() const { return eps_; }
};

/*
 Eq. 12b [2]
 */
class GradWHMatrix : public GradCmplxMolMat
{
protected:
  double kappa_;
  
public:
  GradWHMatrix(int I, int wrt, int ns, int p, double kappa);
  
  void calc_all_vals(shared_ptr<Molecule> mol,
                     shared_ptr<PrecalcBessel> bcalc,
                     shared_ptr<GradHMatrix> dH,
                     shared_ptr<GradFMatrix> dF,
                     shared_ptr<GradLHMatrix> dLH,
                     shared_ptr<GradLHNMatrix> dLHN,
                     shared_ptr<GradLFMatrix> dLF);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  shared_ptr<PrecalcBessel> bcalc,
                  shared_ptr<GradHMatrix> dH,
                  shared_ptr<GradFMatrix> dF,
                  shared_ptr<GradLHMatrix> dLH,
                  shared_ptr<GradLHNMatrix> dLHN,
                  shared_ptr<GradLFMatrix> dLF);
};

/*
 Eq. 11a [2]
 */
class GradFMatrix : public GradCmplxMolMat
{
public:
  GradFMatrix(int I, int wrt, int ns, int p);
  
  void calc_all_vals(shared_ptr<IEMatrix> IE,
                     shared_ptr<GradWFMatrix> dWF);
  
  void calc_val_k(int k, shared_ptr<IEMatrix> IE,
                  shared_ptr<GradWFMatrix> dWF);
  
};

/*
 Eq. 11b [2]
 */
class GradHMatrix : public GradCmplxMolMat
{
protected:
  double kappa_;
  
public:
  GradHMatrix(int I, int wrt,
              int ns, int p, double kappa);
  
  void calc_all_vals(shared_ptr<PrecalcBessel> bcalc,
                     shared_ptr<IEMatrix> IE,
                     shared_ptr<GradWHMatrix> dWH);
  
  void calc_val_k(int k,
                  shared_ptr<PrecalcBessel> bcalc,
                  shared_ptr<IEMatrix> IE,
                  shared_ptr<GradWHMatrix> dWH);
  
};

/*
 Eq. 13 [2]
 */
class GradLFMatrix : public GradCmplxMolMat
{
public:
  GradLFMatrix(int I, int wrt, int ns, int p);
  
  void calc_all_vals(shared_ptr<Molecule> mol,
                     shared_ptr<SHCalc> shcalc,
                     shared_ptr<TMatrix> T,
                     shared_ptr<GradFMatrix> dF,
                     int Mp=-1);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  shared_ptr<SHCalc> shcalc,
                  shared_ptr<TMatrix> T,
                  shared_ptr<GradFMatrix> dF,
                  int Mp=-1);
  
  MyMatrix<Ptx> numeric_reex(int k, int j,
                             shared_ptr<Molecule> mol,
                             shared_ptr<SHCalc> shcalc,
                             shared_ptr<GradFMatrix> dF,
                             int Mp=-1);
  
  // calculate the gradient of h at point P (Eq. S5a)
  Ptx calc_df_P(Pt P, int k, shared_ptr<SHCalc> shcalc,
                shared_ptr<GradFMatrix> dF);
  
};

/*
 Eq. 13 [2]
 */
class GradLHMatrix : public GradCmplxMolMat
{
protected:
  double kappa_;
  
public:
  GradLHMatrix(int I, int wrt, int ns, int p, double kappa);
  
  void calc_all_vals(shared_ptr<Molecule> mol,
                     shared_ptr<PrecalcBessel> bcalc,
                     shared_ptr<SHCalc> shcalc,
                     shared_ptr<TMatrix> T,
                     shared_ptr<GradHMatrix> dH, int Mp=-1);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  shared_ptr<PrecalcBessel> bcalc,
                  shared_ptr<SHCalc> shcalc,
                  shared_ptr<TMatrix> T,
                  shared_ptr<GradHMatrix> dH, int Mp=-1);
  
  MyMatrix<Ptx> numeric_reex(int k, int j,
                             shared_ptr<Molecule> mol,
                             shared_ptr<PrecalcBessel> bcalc,
                             shared_ptr<SHCalc> shcalc,
                             shared_ptr<GradHMatrix> dH,
                             int Mp=-1);
  
  // calculate the gradient of h at point P (Eq. S5a)
  Ptx calc_dh_P(Pt P, int k, shared_ptr<PrecalcBessel> bcalc,
                shared_ptr<SHCalc> shcalc,
                shared_ptr<GradHMatrix> dH);
  
};

/*
 Eq. 13 [2]
 */
class GradLHNMatrix : public GradCmplxMolMat
{
public:
  GradLHNMatrix(int I, int wrt, int ns, int p);
  
  void calc_all_vals(shared_ptr<System> sys, shared_ptr<TMatrix> T,
                     vector<shared_ptr<HMatrix> > H,
                     vector<shared_ptr<GradHMatrix> > dH);
  
  void calc_val_k(int k, shared_ptr<System> sys,
                  shared_ptr<TMatrix> T,
                  vector<shared_ptr<HMatrix> > H,
                  vector<shared_ptr<GradHMatrix> > dH);
  
};

#endif /* Gradsolvmat_h */

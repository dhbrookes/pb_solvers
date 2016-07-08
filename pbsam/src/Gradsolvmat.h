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
  
  cmplx get_mat_knm_d(int k, int n, int m, int d)
  {
    if ( d == 0 )      return mat_[k](n, m+p_).x();
    else if ( d == 1 ) return mat_[k](n, m+p_).y();
    return mat_[k](n, m+p_).z();
  }
  void set_mat_knm(int k, int n, int m, Ptx val)
  { mat_[k].set_val(n, m+p_, val); }
  
  MyMatrix<Ptx> get_mat_k(int k) const { return mat_[k]; }
  void set_mat_k(int k, MyMatrix<Ptx> mat )
  {
    for (int n = 0; n < p_; n++)
      for (int m = -n; m <= n; m++)
        mat_[k].set_val(n, m+p_, mat(n, m+p_));
  }
  
  void add_mat_k(int k, MyMatrix<Ptx> mat )
  {
    for (int n = 0; n < p_; n++)
      for (int m = -n; m <= n; m++)
        mat_[k].set_val(n, m+p_, mat_[k](n, m+p_) + mat(n, m+p_));
  }
  
  void reset_mat();
  
  friend ostream & operator<<(ostream & fout, GradCmplxMolMat & M)
  {
    //    fout << "{{";
    for (int k = 0; k < M.get_ns(); k++)
    {
      fout << "For sphere " << k << endl;
      for (int d = 0; d < 3; d++)
      {
        fout << " Dim: " << d <<  endl;
        for (int n = 0; n < M.get_p(); n++)
        {
          for (int m = 0; m <= n; m++)
          {
            double real = M.get_mat_knm_d( k, n, m, d).real();
            double imag = M.get_mat_knm_d( k, n, m, d).imag();
            if(abs(real) < 1e-15 ) real = 0.0;
            if(abs(imag) < 1e-15 ) imag = 0.0;
            fout << "(" << setprecision(7)<<  real << ", " << imag << ") ";
            //          fout << setprecision(9) << real << ",";
          }
          fout << endl;
        }
        //      fout << "},{" ;
        fout << endl;
      }
      //      fout << "},{" ;
      fout << endl;
    }
    //    fout << "},{" << endl;
    return fout;
  }
  
  void print_kmat(int k)
  {
    cout << "Molecule " << I_ << " For sphere " << k << endl;
    for (int d = 0; d < 3; d++)
    {
      cout << " Dim: " << d <<  endl;
      for (int n = 0; n < get_p(); n++)
      {
        for (int m = 0; m <= n; m++)
        {
          double real = get_mat_knm_d( k, n, m, d).real();
          double imag = get_mat_knm_d( k, n, m, d).imag();
          if(abs(real) < 1e-15 ) real = 0.0;
          if(abs(imag) < 1e-15 ) imag = 0.0;
          cout << setprecision(9) << "(" << real << ", " << imag << ") ";
        }
        cout << endl;
      }
      cout << endl;
    }
    cout << endl;
  }

};

/*
 Base class for gradients of NumericalMatrix objects
 */
class GradNumericalMat
{
protected:
  int p_;
  int I_;
  int wrt_;  // with respect to
  vector<MyMatrix<Ptx> > mat_cmplx_;  // each gradient has 3 parts (use Pt)
  vector<vector<Pt> > mat_;
  
  
  
public:
  GradNumericalMat(int I, int wrt, int ns, int p)
  :p_(p), mat_(ns), mat_cmplx_(ns, MyMatrix<Ptx> (p, 2*p+1)), I_(I), wrt_(wrt)
  {
  }
  
  const int get_wrt() const   { return wrt_; }
  const int get_I() const   { return I_; }
  const int get_p() const   { return p_; }
  const int get_ns() const  { return (int) mat_.size(); }

  Pt get_mat_kh(int k, int h)           { return mat_[k][h]; }
  Ptx get_mat_knm(int k, int n, int m)  { return mat_cmplx_[k](n,m+p_); }
  void set_mat_kh(int k, int h, Pt val) { mat_[k][h] = val; }
  vector<Pt> get_mat_k(int k)           { return mat_[k]; }
  vector<vector<Pt> > get_mat()         { return mat_; }
  int get_mat_k_len(int k)                 { return (int)mat_[k].size(); }
  
  void reset_mat(int k);
  
  friend ostream & operator<<(ostream & fout, GradNumericalMat & M)
  {
    for (int k = 0; k < M.get_ns(); k++)
    {
      fout << "For sphere " << k << endl;
      for (int d = 0; d < 3; d++)
      {
        fout << " Dim: " << d <<  endl;
        for (int h = 0; h < M.get_mat_k_len(k); h++)
        {
          double real;
          if (d == 0)       real = M.get_mat_kh( k, h).x();
          else if (d == 1)  real = M.get_mat_kh( k, h).y();
          else              real = M.get_mat_kh( k, h).z();
          if(abs(real) < 1e-15 ) real = 0.0;
          fout << real << ", ";
        }
        fout << endl;
      }
    }
    return fout;
  }
  
  void print_kmat(int k)
  {
    cout << "Molecule " << I_ << " For sphere " << k << endl;
    for (int d = 0; d < 3; d++)
    {
      cout << " Dim: " << d <<  endl;
      for (int h = 0; h < get_mat_k_len(k); h++)
      {
        double real;
        if (d == 0)       real = get_mat_kh( k, h).x();
        else if (d == 1)  real = get_mat_kh( k, h).y();
        else              real = get_mat_kh( k, h).z();
        if(abs(real) < 1e-15 ) real = 0.0;
        cout << real << ", ";
      }
      cout << endl;
    }
    cout << endl;
  }
  
  void print_analytical(int k)
  {
    cout << "Molecule " << I_ << " For sphere " << k << endl;
    for (int d = 0; d < 3; d++)
    {
      cout << " Dim: " << d <<  endl;
      for (int n = 0; n < get_p(); n++)
      {
        for (int m = 0; m <= n; m++)
        {
          double real, imag;
          if (d == 0)
          {
            real = (get_mat_knm( k, n, m).x()).real();
            imag = (get_mat_knm( k, n, m).x()).imag();
          } else if (d == 1)
          {
            real = (get_mat_knm( k, n, m).y()).real();
            imag = (get_mat_knm( k, n, m).y()).imag();
          } else
          {
            real = (get_mat_knm( k, n, m).z()).real();
            imag = (get_mat_knm( k, n, m).z()).imag();
          }
          
          if(abs(real) < 1e-15 ) real = 0.0;
          if(abs(imag) < 1e-15 ) imag = 0.0;
          cout << setprecision(9)<< "(" << real << ", " << imag << ") ";
        }
        cout << endl;
      }
    }
  }

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
                     shared_ptr<BesselCalc> bcalc,
                     shared_ptr<GradHMatrix> dH,
                     shared_ptr<GradFMatrix> dF,
                     shared_ptr<GradLHMatrix> dLH,
                     shared_ptr<GradLHNMatrix> dLHN,
                     shared_ptr<GradLFMatrix> dLF);
  
  // calcualte values for one sphere in this molecule
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  vector<double> besseli,
                  vector<double> besselk,
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
                     shared_ptr<BesselCalc> bcalc,
                     shared_ptr<GradHMatrix> dH,
                     shared_ptr<GradFMatrix> dF,
                     shared_ptr<GradLHMatrix> dLH,
                     shared_ptr<GradLHNMatrix> dLHN,
                     shared_ptr<GradLFMatrix> dLF);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol,
                  vector<double> besseli,
                  vector<double> besselk,
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
  
  Ptx calc_df_P(Pt P, int k, shared_ptr<SHCalc> shcalc);
  
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
  
  void calc_all_vals(shared_ptr<Molecule> mol,
                     shared_ptr<BesselCalc> bcalc,
                     shared_ptr<IEMatrix> IE,
                     shared_ptr<GradWHMatrix> dWH);
  
  void calc_val_k(int k,
                  vector<double> besseli,
                  shared_ptr<IEMatrix> IE,
                  shared_ptr<GradWHMatrix> dWH);
  
  // calculate the gradient of h at point P (Eq. S5a)
  Ptx calc_dh_P(Pt P, int k, vector<double> besseli,
                shared_ptr<SHCalc> shcalc);
  
  double get_kappa() const { return kappa_; }
  
};

/*
 Eq. 13 [2]
 */
class GradLFMatrix : public GradNumericalMat
{
public:
  GradLFMatrix(int I, int wrt, int ns, int p);
  
  void init(shared_ptr<Molecule> mol, shared_ptr<GradFMatrix> dF,
            shared_ptr<SHCalc> shcalc,
            shared_ptr<ExpansionConstants> _expconst);
  
  void calc_all_vals(shared_ptr<Molecule> mol, vector<int> interpol,
                     shared_ptr<TMatrix> T, shared_ptr<GradFMatrix> dF);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol, vector<int> interpol,
                  shared_ptr<TMatrix> T,
                  shared_ptr<GradFMatrix> dF);
  
//  MyMatrix<Ptx> numeric_reex(int k, int j,
//                             shared_ptr<Molecule> mol,
//                             shared_ptr<SHCalc> shcalc,
//                             shared_ptr<GradFMatrix> dF,
//                             int Mp=-1);
};

/*
 Eq. 13 [2]
 */
class GradLHMatrix : public GradNumericalMat
{
protected:
  double kappa_;
  
public:
  GradLHMatrix(int I, int wrt, int ns, int p, double kappa);
  
  void init(shared_ptr<Molecule> mol, shared_ptr<GradHMatrix> dH,
            shared_ptr<SHCalc> shcalc, shared_ptr<BesselCalc> bcalc,
            shared_ptr<ExpansionConstants> _expconst);
  
  void calc_all_vals(shared_ptr<Molecule> mol, vector<int> interpol,
                     shared_ptr<TMatrix> T, shared_ptr<GradHMatrix> dH);
  
  void calc_val_k(int k, shared_ptr<Molecule> mol, vector<int> interpol,
                  shared_ptr<TMatrix> T, shared_ptr<GradHMatrix> dH);
  
//  MyMatrix<Ptx> numeric_reex(int k, int j,
//                             shared_ptr<Molecule> mol,
//                             vector<double> besseli,
//                             vector<double> besselk,
//                             shared_ptr<SHCalc> shcalc,
//                             shared_ptr<GradHMatrix> dH,
//                             int Mp=-1);
  

  
};

/*
 Eq. 13 [2]
 */
class GradLHNMatrix : public GradCmplxMolMat
{
public:
  GradLHNMatrix(int I, int wrt, int ns, int p);
  
  void calc_all_vals(shared_ptr<System> sys, shared_ptr<TMatrix> T,
                     vector<shared_ptr<GradCmplxMolMat> > gradT_A,
                     vector<shared_ptr<GradHMatrix> > dH);
  
  void calc_val_k(int k, shared_ptr<System> sys,
                  shared_ptr<TMatrix> T,
                  vector<shared_ptr<GradCmplxMolMat> > gradT_A,
                  vector<shared_ptr<GradHMatrix> > dH);
  
};

#endif /* Gradsolvmat_h */

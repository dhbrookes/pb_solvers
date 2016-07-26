//
//  Gradsolvmat.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Gradsolvmat.h"


void GradCmplxMolMat::reset_mat()
{
  for (int k = 0; k < mat_.size(); k++)
  {
    for (int i = 0; i < mat_[k].get_nrows(); i++)
    {
      for (int j = 0; j < mat_[k].get_ncols(); j++)
      {
        mat_[k].set_val(i, j, Ptx());
      }
    }
  }
}


void GradNumericalMat::reset_mat(int k)
{
  for (int i = 0; i < mat_cmplx_[k].get_nrows(); i++)
  {
    for (int j = 0; j < mat_cmplx_[k].get_ncols(); j++)
    {
      mat_cmplx_[k].set_val(i, j, cmplx(0, 0));
    }
  }
  
}

GradWFMatrix::GradWFMatrix(int I, int wrt, int ns, int p,
                           double eps_in, double eps_out,
                           double kappa)
:GradCmplxMolMat(I, wrt, ns, p), eps_(eps_in/eps_out), kappa_(kappa)
{
}



void GradWFMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                                 shared_ptr<BesselCalc> bcalc,
                                 shared_ptr<GradHMatrix> dH,
                                 shared_ptr<GradFMatrix> dF,
                                 shared_ptr<GradLHMatrix> dLH,
                                 shared_ptr<GradLHNMatrix> dLHN,
                                 shared_ptr<GradLFMatrix> dLF)
{
  vector<double> besseli, besselk;
  for (int k = 0; k < get_ns(); k++)
  {
    besseli = bcalc->calc_mbfI(p_+1, kappa_*mol->get_ak(k));
    besselk = bcalc->calc_mbfK(p_+1, kappa_*mol->get_ak(k));
    calc_val_k(k, mol, besseli, besselk, dH, dF, dLH, dLHN, dLF);
  }
}

void GradWFMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              vector<double> besseli,
                              vector<double> besselk,
                              shared_ptr<GradHMatrix> dH,
                              shared_ptr<GradFMatrix> dF,
                              shared_ptr<GradLHMatrix> dLH,
                              shared_ptr<GradLHNMatrix> dLHN,
                              shared_ptr<GradLFMatrix> dLF)
{
  double ak = mol->get_ak(k);
  double exp_ka = exp(-kappa_*ak);
  double bin, bin1, bkn, bkn1;
  Ptx tot, val1, val2, val3, val4;
  double inner;
  for (int n = 0; n < p_; n++)
  {
    bin = besseli[n];
    bin1 = besseli[n+1];
    bkn = besselk[n];
    bkn1 = besselk[n+1];
    for (int m = -n; m < n+1; m++)
    {
      val1=dH->get_mat_knm(k,n,m)*exp_ka*(n*bkn-(2*n+1)*bkn1);
      val2=dF->get_mat_knm(k,n,m)*(2*n+1-n*eps_);
      val3=dLF->get_mat_knm(k, n, m)*n*eps_*ak;
      inner=(n*bin)+((bin1*pow(kappa_*ak,2))/(double)(2*n+3));
      val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m)) * inner * ak;
      set_mat_knm(k, n, m, val1+val2-val3+val4);
    }
  }
}


GradWHMatrix::GradWHMatrix(int I, int wrt,
                           int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}


void GradWHMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                             shared_ptr<BesselCalc> bcalc,
                             shared_ptr<GradHMatrix> dH,
                             shared_ptr<GradFMatrix> dF,
                             shared_ptr<GradLHMatrix> dLH,
                             shared_ptr<GradLHNMatrix> dLHN,
                             shared_ptr<GradLFMatrix> dLF)
{
  vector<double> besseli, besselk;
  for (int k = 0; k < get_ns(); k++)
  {
    besseli = bcalc->calc_mbfI(p_+1, kappa_*mol->get_ak(k));
    besselk = bcalc->calc_mbfK(p_+1, kappa_*mol->get_ak(k));
    calc_val_k(k, mol, besseli, besselk, dH, dF, dLH, dLHN, dLF);
  }
}


void GradWHMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              vector<double> besseli,
                              vector<double> besselk,
                              shared_ptr<GradHMatrix> dH,
                              shared_ptr<GradFMatrix> dF,
                              shared_ptr<GradLHMatrix> dLH,
                              shared_ptr<GradLHNMatrix> dLHN,
                              shared_ptr<GradLFMatrix> dLF)
{
  double ak = mol->get_ak(k);
  double exp_ka = exp(-kappa_*ak);
  double bin, bkn, inner;
  Ptx tot, val1, val2, val3, val4;
  
//  cout << "Inside WH" << endl;
//  cout << "This is dLH ";
//  dLH->print_analytical(k);
//  cout << "This is dLHN ";
//  dLHN->print_kmat(k);
  
//  cout << "This is dLF ";
//  dLF->print_analytical(k);
//  cout << "This is dF ";
//  dF->print_kmat(k);
  
  for (int n = 0; n < p_; n++)
  {
    bin = besseli[n];
    bkn = besselk[n];
    for (int m = -n; m < n+1; m++)
    {
      inner = ((2.0*n+1.0) / bin) - (exp_ka * bkn);
      val1 = dH->get_mat_knm(k, n, m) * inner;
      val2 = dF->get_mat_knm(k, n, m);
      val3 = dLF->get_mat_knm(k, n, m) * ak;
      val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m))*ak*bin;
      tot = val1+val2+val3-val4;
      set_mat_knm(k, n, m, tot);
    }
  }
}

GradFMatrix::GradFMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradFMatrix::calc_all_vals(shared_ptr<IEMatrix> IE,
                            shared_ptr<GradWFMatrix> dWF)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, IE, dWF);
}

void GradFMatrix::calc_val_k(int k, shared_ptr<IEMatrix> IE,
                             shared_ptr<GradWFMatrix> dWF)
{
  int ct;
  vector<MyMatrix<double> > df_in(3, MyMatrix<double> (p_*p_, 1));
  vector<MyMatrix<double> > df_out(3, MyMatrix<double> (p_*p_, 1));
  for (int d = 0; d < 3; d++)  // 3 dimensions
  {
    ct = 0;
    for (int l = 0; l < p_; l++)  // rows in old matrix
    {
      for (int s = 0; s < l+1; s++)  //columns in old matrix
      {
        df_in[d](ct,0) = dWF->get_mat_knm_d(k, l, s, d).real();
        ct++;
        if ( s > 0 )
        {
          df_in[d](ct,0) = dWF->get_mat_knm_d(k, l, s, d).imag();
          ct++;
        }
      }
    }
    //TODO: replace with matMul
    df_out[d] = IE->get_IE_k(k) * df_in[d];  // Matrix vector multiplication
  }

  for (int d = 0; d < 3; d++)  // 3 dimensions
  {
    int ctr(0);
    for (int n = 0; n < p_; n++)  // rows in new matrix
    {
      double scl = 1.0 / (double)(2.0*n+1.0);
      for (int m = 0; m < n+1; m++)  // columns in new matrix
      {
        double fRe, fIm;
        fRe = df_out[d](ctr,0);
        ctr++;
        
        if ( m > 0 )
        {
          fIm = df_out[d](ctr,0);
          ctr++;
        } else
          fIm = 0.0;
        
        set_mat_knm_d(k,n,m,d, scl*complex<double> (fRe,fIm));
        if ( m > 0 ) set_mat_knm_d(k,n,-m,d,scl*complex<double> (fRe,-fIm));
      }
    }
  }
}

Ptx GradFMatrix::calc_df_P(Pt P, int k, shared_ptr<SHCalc> shcalc)
{
  shcalc->calc_sh(P.theta(), P.phi());
  Ptx df = Ptx();
  Ptx in;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = get_mat_knm(k,n,m)*shcalc->get_result(n,m)*((2*n+1)/(4*M_PI));
      df += in;
    }
  }
  return df;
}


GradHMatrix::GradHMatrix(int I, int wrt, int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}

void GradHMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                                shared_ptr<BesselCalc> bcalc,
                                shared_ptr<IEMatrix> IE,
                                shared_ptr<GradWHMatrix> dWH)
{
  vector<double> besseli;
  for (int k = 0; k < get_ns(); k++)
  {
    besseli = bcalc->calc_mbfI(p_, kappa_*mol->get_ak(k));
    calc_val_k(k, besseli, IE, dWH);
  }
}

void GradHMatrix::calc_val_k(int k, vector<double> besseli,
                             shared_ptr<IEMatrix> IE,
                             shared_ptr<GradWHMatrix> dWH)
{
  int ct, d, l, s, n, m;
  MyMatrix<double> dh_in(p_*p_, 1);
  vector<MyMatrix<double> > dh_out(3, MyMatrix<double> (p_*p_, 1));
  for (d = 0; d < 3; d++)  // 3 dimensions
  {
    ct = 0;
    for (l = 0; l < p_; l++)  // rows in old matrix
    {
      for (s = 0; s < l+1; s++)  //columns in old matrix
      {
        dh_in(ct,0) = dWH->get_mat_knm_d(k, l, s, d).real();
        ct++;
        if ( s > 0 )
        {
          dh_in(ct,0) = dWH->get_mat_knm_d(k, l, s, d).imag();
          ct++;
        }
      }
    }
    //TODO: replace with matMul
    dh_out[d] = IE->get_IE_k(k) * dh_in;  // Matrix vector multiplication
  }
  
//  
//  cout << "~~~~~~ output ~~~~~~~~~~~~~~" << endl;
//  for (int d = 0; d < 3; d++)
//  {
//    ct = 0;
//    cout << " Dim: " << d <<  endl;
//    for (int n = 0; n < p_; n++)
//    {
//      for (int m = 0; m < p_; m++)
//      {
//        double real = dh_out[d](ct, 0); //.real();
//        if(abs(real) < 1e-15 ) real = 0.0;
//        cout << setprecision(7)<<  real << ", ";
//        ct++;
//      }
//      cout << endl;
//    }
//    cout << endl;
//  }
//  
  
  for (d = 0; d < 3; d++)  // 3 dimensions
  {
    int ctr(0);
    for (n = 0; n < p_; n++)  // rows in new matrix
    {
      double scl = besseli[n] / (double)(2.0*n+1.0);
      for (m = 0; m < n+1; m++)  // columns in new matrix
      {
        double dhR, dhI;
        dhR = dh_out[d](ctr,0);
        ctr++;
        
        if ( m > 0 )
        {
          dhI = dh_out[d](ctr,0);
          ctr++;
        } else
          dhI = 0.0;
        
        set_mat_knm_d(k,n,m,d,scl*complex<double> (dhR,dhI));
        if ( m > 0 ) set_mat_knm_d(k,n,-m,d,scl*complex<double> (dhR,-dhI));
      }
    }
  }
}

GradLFMatrix::GradLFMatrix(int I, int wrt, int ns, int p)
:GradNumericalMat(I, wrt, ns, p)
{
}

void GradLFMatrix::init_k(int k, shared_ptr<Molecule> mol,
                          shared_ptr<GradFMatrix> dF,
                          shared_ptr<SHCalc> shcalc,
                          shared_ptr<ExpansionConstants> _expconst)
{
  double rl, im;
  
  for (int k=0; k< mol->get_ns(); k++)
  {
    vector <int> exp_pts = mol->get_gdpt_expj(k);
    mat_[k].resize( (int) exp_pts.size() );
    double dA = 4 * M_PI / (double) mol->get_gridj(k).size();
    
    for (int h = 0; h < exp_pts.size(); h++)
    {
      vector<double> val(3, 0.0);
      Pt q = mol->get_gridjh(k, exp_pts[h]);
      shcalc->calc_sh(q.theta(), q.phi());
      
      for (int d = 0; d < 3; d++)
      {
        for (int n=0; n<p_; n++)
        {
          for (int m = -n; m <= n; m++)
          {
            rl = (shcalc->get_result(n,m).real()*
                  dF->get_mat_knm_d(k,n,m,d).real());
            im = (shcalc->get_result(n,m).imag()*
                  dF->get_mat_knm_d(k,n,m,d).imag());
            val[d] += dA*_expconst->get_const1_l(n) * ( rl + im );
          }
        }
      }
      set_mat_kh(k, h, Pt(val[0], val[1], val[2]));
    }
  }
}

void GradLFMatrix::calc_all_vals(shared_ptr<Molecule> mol, vector<int> interpol,
                                 shared_ptr<TMatrix> T,
                                 shared_ptr<GradFMatrix> dF)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, interpol, T, dF);
}

void GradLFMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              vector<int> interpol,
                              shared_ptr<TMatrix> T,
                              shared_ptr<GradFMatrix> dF)
{
  MyMatrix<Ptx> reex, inner(p_, 2*p_+1);
  
  for (int j = 0; j < get_ns(); j++)
  {
    if (j == k) continue;
    if ((interpol[k] != 0) || (interpol[j] != 0)) continue;
    
    if (T->is_analytic(I_, k, I_, j))
    {
//      cout << "LF Reexp Analyt " << I_ << " sph " << j << " to " <<
//      I_ << " sph " << k << endl;
//      dF->print_kmat(j);
      reex = T->re_expand_gradX(dF->get_mat_k(j), I_, k, I_, j, true);
      
//      cout << "~~~~~~ output ~~~~~~~~~~~~~~" << endl;
//      for (int d = 0; d < 3; d++)
//      {
//        cout << " Dim: " << d <<  endl;
//        for (int n = 0; n < p_; n++)
//        {
//          for (int m = 0; m <= n; m++)
//          {
//            double real = reex( n, m+p_).get_cart(d).real();
//            double imag = reex( n, m+p_).get_cart(d).imag();
//            if(abs(real) < 1e-15 ) real = 0.0;
//            if(abs(imag) < 1e-15 ) imag = 0.0;
//            cout << "(" << setprecision(7)<<  real << ", " << imag << ") ";
//          }
//          cout << endl;
//        }
//        cout << endl;
//      }

    }
    else
    {
      // From (I,j) -> (I,k) is re_exp_numeric( , I, k, I, j)
//      cout << "LF Reexp Numer " << I_ << " sph " << j << " to " <<
//      I_ << " sph " << k << endl;
//      print_kmat(j);
      reex = T->re_expandgradX_numeric(get_mat(), I_, k, I_, j, 0.0 );
    }
    
    inner += reex;
  }
  mat_cmplx_[k] = inner;
}


GradLHMatrix::GradLHMatrix(int I, int wrt, int ns, int p, double kappa)
:GradNumericalMat(I, wrt, ns, p), kappa_(kappa)
{
}

void GradLHMatrix::init_k(int k, shared_ptr<Molecule> mol,
                          shared_ptr<GradHMatrix> dH,
                          shared_ptr<SHCalc> shcalc,
                          shared_ptr<BesselCalc> bcalc,
                          shared_ptr<ExpansionConstants> _expconst)
{
  double rl, im;

  vector <int> exp_pts = mol->get_gdpt_expj(k);
  mat_[k].resize( (int) exp_pts.size() );
  double dA = 4 * M_PI / (double) mol->get_gridj(k).size();
  
  for (int h = 0; h < exp_pts.size(); h++)
  {
    vector<double> val(3, 0.0);
    Pt q = mol->get_gridjh(k, exp_pts[h]);
    shcalc->calc_sh(q.theta(), q.phi());
    vector<double> bessI = bcalc->calc_mbfI(p_+1, kappa_*q.r());
    
    for (int d = 0; d < 3; d++)
    {
      for (int n=0; n<p_; n++)
      {
        for (int m = -n; m <= n; m++)
        {
          rl = (shcalc->get_result(n,m).real()*
                dH->get_mat_knm_d(k,n,m,d).real());
          im = (shcalc->get_result(n,m).imag()*
                dH->get_mat_knm_d(k,n,m,d).imag());
          val[d] += dA*_expconst->get_const1_l(n) * ( rl + im )/bessI[n];
        }
      }
    }
    set_mat_kh(k, h, Pt(val[0], val[1], val[2]));
  }
}

void GradLHMatrix::calc_all_vals(shared_ptr<Molecule> mol, vector<int> interpol,
                                 shared_ptr<TMatrix> T,
                                 shared_ptr<GradHMatrix> dH)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, interpol, T, dH);
}

void GradLHMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              vector<int> interpol,
                              shared_ptr<TMatrix> T,
                              shared_ptr<GradHMatrix> dH)
{
  MyMatrix<Ptx> reex, inner(p_, 2*p_+1);

  for (int j = 0; j < get_ns(); j++)
  {
    if (j == k) continue;
    if ((interpol[k] != 0) || (interpol[j] != 0)) continue;
    
    if (T->is_analytic(I_, k, I_, j))
    {
//      cout << "LH Reexp Analyt " << I_ << " sph " << j << " to " <<
//      I_ << " sph " << k << endl;
//      dH->print_kmat(j);
      reex = T->re_expand_gradX(dH->get_mat_k(j), I_, k, I_, j);
    }
    else
    {
//      cout << "LH Reexp Numer " << I_ << " sph " << j << " to " <<
//      I_ << " sph " << k << endl;
//      print_kmat(j);
      reex = T->re_expandgradX_numeric(get_mat(), I_, k, I_, j, kappa_);
    }
    
//    cout << "~~~~~~ output ~~~~~~~~~~~~~~" << endl;
//    for (int d = 0; d < 3; d++)
//    {
//      cout << " Dim: " << d <<  endl;
//      for (int n = 0; n < p_; n++)
//      {
//        for (int m = 0; m <= n; m++)
//        {
//          double real = reex( n, m+p_).get_cart(d).real();
//          double imag = reex( n, m+p_).get_cart(d).imag();
//          if(abs(real) < 1e-15 ) real = 0.0;
//          if(abs(imag) < 1e-15 ) imag = 0.0;
//          cout << "(" << setprecision(7)<<  real << ", " << imag << ") ";
//        }
//        cout << endl;
//      }
//      cout << endl;
//    }
    inner += reex;
  }
  mat_cmplx_[k] = inner;
}


Ptx GradHMatrix::calc_dh_P(Pt P, int k,
                            vector<double> besseli,
                            shared_ptr<SHCalc> shcalc)
{
  shcalc->calc_sh(P.theta(), P.phi());
  Ptx dh = Ptx();
  Ptx in;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = get_mat_knm(k,n,m)*shcalc->get_result(n,m)*((2*n+1)/(4*M_PI));
      in *= 1/besseli[n];
      dh += in;
    }
  }
  return dh;
}


GradLHNMatrix::GradLHNMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradLHNMatrix::calc_all_vals(shared_ptr<System> sys,
                                  shared_ptr<TMatrix> T,
                                  vector<shared_ptr<GradCmplxMolMat> > gradT_A,
                                  vector<shared_ptr<GradHMatrix> > dH)
{
  MyMatrix<Ptx> lhn_k, inner;
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, sys, T, gradT_A, dH);
}

void GradLHNMatrix::calc_val_k(int k, shared_ptr<System> sys,
                               shared_ptr<TMatrix> T,
                               vector<shared_ptr<GradCmplxMolMat> > gradT_A,
                               vector<shared_ptr<GradHMatrix> > dH,
                               double dist_cut)
{
  Pt Ik, Mm;
  double aIk, aMm, interPolcut = dist_cut;
  MyMatrix<Ptx> lhn_k(p_, 2*p_+1), inner;
  lhn_k = gradT_A[I_]->get_mat_k(k);
  
  if ( fabs(interPolcut-100) < 1e-2)
  {
  cout << "~~~~~~ " << I_ << " and k " << k <<" gT * A ~~~~~~~~~~~~~~" << endl;
  for (int d = 0; d < 3; d++)
  {
    cout << " Dim: " << d <<  endl;
    for (int n = 0; n < p_; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double real = lhn_k( n, m+p_).get_cart(d).real();
        double imag = lhn_k( n, m+p_).get_cart(d).imag();
        if(abs(real) < 1e-15 ) real = 0.0;
        if(abs(imag) < 1e-15 ) imag = 0.0;
        cout << "(" << setprecision(7)<<  real << ", " << imag << ") ";
      }
      cout << endl;
    }
  }
  cout << "~~~~~~~~~~~~~~~~~~~~" << endl;
  }
  
  Ik = sys->get_centerik(I_, k);
  aIk = sys->get_aik(I_, k);
  
  for (int M=0; M < sys->get_n(); M++)
  {
    if (M == I_) continue;
    for (int m=0; m < sys->get_Ns_i(M); m++)
    {
      Mm = sys->get_centerik(M, m);
      aMm = sys->get_aik(M, m);
      
      if ( sys->get_pbc_dist_vec_base(Ik, Mm).norm() < (interPolcut+aIk+aMm))
      {
        if ( fabs(interPolcut-100) < 1e-2)
        {
        cout << "Reex " << I_ << " and sph " << k  << " with xforms from "
        << M << " sph " << m << endl;
        dH[M]->print_kmat(m);
        }
        inner = T->re_expand_gradX(dH[M]->get_mat_k(m), I_, k, M, m);
        
        if ( fabs(interPolcut-100) < 1e-2)
        {
        cout << "~~~~~~ output ~~~~~~~~~~~~~~" << endl;
        for (int d = 0; d < 3; d++)
        {
          cout << " Dim: " << d <<  endl;
          for (int n = 0; n < p_; n++)
          {
            for (int m = 0; m <= n; m++)
            {
              double real = inner( n, m+p_).get_cart(d).real();
              double imag = inner( n, m+p_).get_cart(d).imag();
              if(abs(real) < 1e-15 ) real = 0.0;
              if(abs(imag) < 1e-15 ) imag = 0.0;
              cout << "(" << setprecision(7)<<  real << ", " << imag << ") ";
            }
            cout << endl;
          }
          cout << endl;
        }
        }
        
        lhn_k += inner;
      }
    }
    mat_[k] = lhn_k;
  }
}

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


GradWFMatrix::GradWFMatrix(int I, int wrt, int ns, int p,
                           double eps_in, double eps_out,
                           double kappa)
:GradCmplxMolMat(I, wrt, ns, p), eps_(eps_in/eps_out), kappa_(kappa)
{
}


void GradWFMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                                 shared_ptr<PrecalcBessel> bcalc,
                                 shared_ptr<GradHMatrix> dH,
                                 shared_ptr<GradFMatrix> dF,
                                 shared_ptr<GradLHMatrix> dLH,
                                 shared_ptr<GradLHNMatrix> dLHN,
                                 shared_ptr<GradLFMatrix> dLF)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, bcalc, dH, dF, dLH, dLHN, dLF);
}

void GradWFMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              shared_ptr<PrecalcBessel> bcalc,
                              shared_ptr<GradHMatrix> dH,
                              shared_ptr<GradFMatrix> dF,
                              shared_ptr<GradLHMatrix> dLH,
                              shared_ptr<GradLHNMatrix> dLHN,
                              shared_ptr<GradLFMatrix> dLF)
{
  double ak = mol->get_ak(k);
  double exp_ka = exp(-kappa_*ak);
  double bin, bin1, bkn, bkn1;
  Ptx val1, val2, val3, val4;
  double inner;
  for (int n = 0; n < p_; n++)
  {
    bin = bcalc->get_in_Ik(I_, k, n);
    bin1 = bcalc->get_in_Ik(I_, k, n+1);
    bkn = bcalc->get_kn_Ik(I_, k, n);
    bkn1 = bcalc->get_kn_Ik(I_, k, n+1);
    
    for (int m = -n; m < n+1; m++)
    {
      val1=dH->get_mat_knm(k,n,m)*exp_ka*(n*bkn-(2*n+1)*bkn1);
      val2=dF->get_mat_knm(k,n,m)*(2*n+1-n*eps_);
      val3=dLF->get_mat_knm(k, n, m)*n*eps_*ak;
      inner=(n*bin)+((bin1*kappa_*ak)/(2*n+3));
      inner+=(bin1*kappa_*ak)/(2*n+3);
      val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m)) * inner;
      set_mat_knm(k, n, m, val1+val2+val3+val4);
    }
  }
}


GradWHMatrix::GradWHMatrix(int I, int wrt,
                           int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}


void GradWHMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                             shared_ptr<PrecalcBessel> bcalc,
                             shared_ptr<GradHMatrix> dH,
                             shared_ptr<GradFMatrix> dF,
                             shared_ptr<GradLHMatrix> dLH,
                             shared_ptr<GradLHNMatrix> dLHN,
                             shared_ptr<GradLFMatrix> dLF)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, bcalc, dH, dF, dLH, dLHN, dLF);
}


void GradWHMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              shared_ptr<PrecalcBessel> bcalc,
                              shared_ptr<GradHMatrix> dH,
                              shared_ptr<GradFMatrix> dF,
                              shared_ptr<GradLHMatrix> dLH,
                              shared_ptr<GradLHNMatrix> dLHN,
                              shared_ptr<GradLFMatrix> dLF)
{
  double ak = mol->get_ak(k);
  double exp_ka = exp(-kappa_*ak);
  double bin, bkn, inner;
  Ptx val1, val2, val3, val4;
  for (int n = 0; n < p_; n++)
  {
    bin = bcalc->get_in_Ik(I_, k, n);
    bkn = bcalc->get_kn_Ik(I_, k, n);
    for (int m = -n; m < n+1; m++)
    {
      inner = ((2 * n + 1) / bin) - (exp_ka * bkn);
      val1 = dH->get_mat_knm(k, n, m) * inner;
      val2 = dF->get_mat_knm(k, n, m);
      val3 = dLF->get_mat_knm(k, n, m) * ak;
      val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m))*ak*bin;
      set_mat_knm(k, n, m, val1+val2+val3+val4);
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
  Ptx in, in2;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = Ptx();
      //compute inner product
      for (int l = 0; l < p_; l++)
      {
        for (int s = 0; s < p_; s++)
        {
          in2=dWF->get_mat_knm(k,l,s);
          //*IE->get_IE_k_nm_ls(k,n,m,l,s).real();
          in += in;
        }
      }
      set_mat_knm(k, n, m, in);
    }
  }
}

GradHMatrix::GradHMatrix(int I, int wrt, int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}

void GradHMatrix::calc_all_vals(shared_ptr<PrecalcBessel> bcalc,
                                shared_ptr<IEMatrix> IE,
                                shared_ptr<GradWHMatrix> dWH)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, bcalc, IE, dWH);
}

void GradHMatrix::calc_val_k(int k, shared_ptr<PrecalcBessel> bcalc,
                             shared_ptr<IEMatrix> IE,
                             shared_ptr<GradWHMatrix> dWH)
{
  Ptx in, in2;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = Ptx();
      //compute inner product
      for (int l = 0; l < p_; l++)
      {
        for (int s = 0; s < p_; s++)
        {
          in2 = dWH->get_mat_knm(k,l,s);
          //*IE->get_IE_k_nm_ls(k,n,m,l,s);
          in += in2;
        }
      }
      set_mat_knm(k, n, m, in * bcalc->get_in_Ik(I_, k, n));
    }
  }
}

GradLFMatrix::GradLFMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradLFMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                                 shared_ptr<SHCalc> shcalc,
                                 shared_ptr<TMatrix> T,
                                 shared_ptr<GradFMatrix> dF, int Mp)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, shcalc, T, dF, Mp);
}

void GradLFMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              shared_ptr<SHCalc> shcalc,
                              shared_ptr<TMatrix> T,
                              shared_ptr<GradFMatrix> dF,
                              int Mp)
{
  MyMatrix<Ptx> reex;
  MyMatrix<Ptx> inner = MyMatrix<Ptx> (p_, 2*p_+1);
  for (int j = 0; j < get_ns(); j++)
  {
    if (j == k) continue;
    if (T->is_analytic(I_, k, I_, j))
    {
      reex = T->re_expand_gradX(dF->get_mat_k(k), I_, k, I_, j);
    }
    else
    {
      reex = numeric_reex(k, j, mol, shcalc, dF);
    }
    inner += reex;
  }
  mat_[k] = inner;
}

MyMatrix<Ptx> GradLFMatrix::numeric_reex(int k, int j,
                                         shared_ptr<Molecule> mol,
                                         shared_ptr<SHCalc> shcalc,
                                         shared_ptr<GradFMatrix> dF,
                                         int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, mol->get_ak(j));
  Ptx fpj, inner;
  Pt rb_k;
  MyMatrix<Ptx> dLF_Ik (p_, 2 * p_ + 1);
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      for (int b = 0; b < Mp; b++)
      {
        fpj = calc_df_P(sph_grid[b], k, shcalc, dF);
        fpj *= (4 * M_PI) / sph_grid.size();
        
        rb_k = sph_grid[b]-(mol->get_centerk(k)-mol->get_centerk(j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        
        inner = fpj * (1/rb_k.r());
        inner *= pow(mol->get_ak(k)/ rb_k.r(), n);
        inner *= conj(shcalc->get_result(n, m));
        dLF_Ik.set_val(n, m, inner);
      }
    }
  }
  return dLF_Ik;
}

Ptx GradLFMatrix::calc_df_P(Pt P, int k, shared_ptr<SHCalc> shcalc,
                            shared_ptr<GradFMatrix> dF)
{
  shcalc->calc_sh(P.theta(), P.phi());
  Ptx df = Ptx();
  Ptx in;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = dF->get_mat_knm(k,n,m)*shcalc->get_result(n,m)*((2*n+1)/(4*M_PI));
      df += in;
    }
  }
  return df;
}


GradLHMatrix::GradLHMatrix(int I, int wrt, int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}

void GradLHMatrix::calc_all_vals(shared_ptr<Molecule> mol,
                                 shared_ptr<PrecalcBessel> bcalc,
                                 shared_ptr<SHCalc> shcalc,
                                 shared_ptr<TMatrix> T,
                                 shared_ptr<GradHMatrix> dH, int Mp)
{
  for (int k = 0; k < get_ns(); k++)
    calc_val_k(k, mol, bcalc, shcalc, T, dH);
}

void GradLHMatrix::calc_val_k(int k, shared_ptr<Molecule> mol,
                              shared_ptr<PrecalcBessel> bcalc,
                              shared_ptr<SHCalc> shcalc,
                              shared_ptr<TMatrix> T,
                              shared_ptr<GradHMatrix> dH, int Mp)
{
  MyMatrix<Ptx> reex, inner;
  inner = MyMatrix<Ptx> (p_, 2*p_+1);
  for (int j = 0; j < get_ns(); j++)
  {
    if (j == k) continue;
    if (T->is_analytic(I_, k, I_, j))
    {
      reex = T->re_expand_gradX(dH->get_mat_k(k), I_, k, I_, j);
    }
    else
    {
      reex = numeric_reex(k, j, mol, bcalc, shcalc, dH, Mp);
    }
    inner += reex;
  }
  mat_[k] = inner;
}

MyMatrix<Ptx> GradLHMatrix::numeric_reex(int k, int j,
                                         shared_ptr<Molecule> mol,
                                         shared_ptr<PrecalcBessel> bcalc,
                                         shared_ptr<SHCalc> shcalc,
                                         shared_ptr<GradHMatrix> dH,
                                         int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, mol->get_ak(j));
  Ptx hpj, inner;
  Pt rb_k;
  MyMatrix<Ptx> dLH_Ik (p_, 2 * p_ + 1);
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      for (int b = 0; b < Mp; b++)
      {
        hpj = calc_dh_P(sph_grid[b], k, bcalc, shcalc, dH);
        hpj *= (4 * M_PI) / sph_grid.size();
        
        rb_k = sph_grid[b]-(mol->get_centerk(k)-mol->get_centerk(j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        
        inner = hpj * (1/rb_k.r());
        inner *= pow(mol->get_ak(k)/ rb_k.r(), n);
        inner *= conj(shcalc->get_result(n, m));
        inner *= exp(-kappa_*rb_k.r());
        inner *= bcalc->get_kn_Ik(I_, k, n);
        
        dLH_Ik.set_val(n, m, inner);
      }
    }
  }
  return dLH_Ik;
}

Ptx GradLHMatrix::calc_dh_P(Pt P, int k,
                            shared_ptr<PrecalcBessel> bcalc,
                            shared_ptr<SHCalc> shcalc,
                            shared_ptr<GradHMatrix> dH)
{
  shcalc->calc_sh(P.theta(), P.phi());
  Ptx dh = Ptx();
  Ptx in;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      in = dH->get_mat_knm(k,n,m)*shcalc->get_result(n,m)*((2*n+1)/(4*M_PI));
      in *= 1/bcalc->get_in_Ik(I_, k, n);
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
                                  vector<shared_ptr<HMatrix> > H,
                                  vector<shared_ptr<GradHMatrix> > dH)
{
  MyMatrix<Ptx> lhn_k, inner;
  for (int k = 0; k < get_ns(); k++)
  {
    lhn_k = MyMatrix<Ptx> (p_, 2*p_+1);
    for (int M=0; M < sys->get_n(); M++)
    {
      for (int m=0; m < sys->get_Ns_i(M); m++)
      {
        inner = T->re_expandX_gradT(H[M]->get_mat_k(m), I_, k, M, m);
        lhn_k += inner;
        
        inner = T->re_expand_gradX(dH[M]->get_mat_k(m), I_, k, M, m);
        lhn_k += inner;
      }
      mat_[k] = lhn_k;
    }
  }
}

void GradLHNMatrix::calc_val_k(int k, shared_ptr<System> sys,
                               shared_ptr<TMatrix> T,
                               vector<shared_ptr<HMatrix> > H,
                               vector<shared_ptr<GradHMatrix> > dH)
{
  MyMatrix<Ptx> lhn_k, inner;
  lhn_k = MyMatrix<Ptx> (p_, 2*p_+1);
  for (int M=0; M < sys->get_n(); M++)
  {
    for (int m=0; m < sys->get_Ns_i(M); m++)
    {
      inner = T->re_expandX_gradT(H[M]->get_mat_k(m), I_, k, M, m);
      lhn_k += inner;
      
      inner = T->re_expand_gradX(dH[M]->get_mat_k(m), I_, k, M, m);
      lhn_k += inner;
    }
    mat_[k] = lhn_k;
  }
}

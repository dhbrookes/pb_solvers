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


void GradWFMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                             shared_ptr<GradHMatrix> dH,
                             shared_ptr<GradFMatrix> dF,
                             shared_ptr<GradLHMatrix> dLH,
                             shared_ptr<GradLHNMatrix> dLHN,
                             shared_ptr<GradLFMatrix> dLF)
{
  double ak;
  double exp_ka;
  Ptx val1, val2, val3, val4;
  double inner;
  vector<double> besselk, besseli;
  for (int k = 0; k < get_ns(); k++)
  {
    ak = mol.get_ak(k);
    besselk = bcalc->calc_mbfK(p_+1, kappa_*ak);
    besseli = bcalc->calc_mbfI(p_+1, kappa_*ak);
    exp_ka = exp(-kappa_*ak);
    for (int n = 0; n < p_; n++)
    {
      for (int m = -n; m < n+1; m++)
      {
        val1=dH->get_mat_knm(k,n,m)*exp_ka*(n*besselk[n]-(2*n+1)*besselk[n+1]);
        val2=dF->get_mat_knm(k,n,m)*(2*n+1-n*eps_);
        val3=dLF->get_mat_knm(k, n, m)*n*eps_*ak;
        inner=(n*besseli[n])+((besseli[n+1]*kappa_*ak)/(2*n+3));
        inner+=(besseli[n+1]*kappa_*ak)/(2*n+3);
        val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m)) * inner;
        set_mat_knm(k, n, m, val1+val2+val3+val4);
      }
    }
  }
}


GradWHMatrix::GradWHMatrix(int I, int wrt, int ns, int p, double kappa)
:GradCmplxMolMat(I, wrt, ns, p), kappa_(kappa)
{
}


void GradWHMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                             shared_ptr<GradHMatrix> dH,
                             shared_ptr<GradFMatrix> dF,
                             shared_ptr<GradLHMatrix> dLH,
                             shared_ptr<GradLHNMatrix> dLHN,
                             shared_ptr<GradLFMatrix> dLF)
{
  double ak;
  double exp_ka;
  Ptx val1, val2, val3, val4;
  double inner;
  vector<double> besselk, besseli;
  for (int k = 0; k < get_ns(); k++)
  {
    ak = mol.get_ak(k);
    besselk = bcalc->calc_mbfK(p_+1, kappa_*ak);
    besseli = bcalc->calc_mbfI(p_+1, kappa_*ak);
    exp_ka = exp(-kappa_*ak);
    for (int n = 0; n < p_; n++)
    {
      for (int m = -n; m < n+1; m++)
      {
        inner = ((2 * n + 1) / besseli[n]) - (exp_ka * besselk[n]);
        val1 = dH->get_mat_knm(k, n, m) * inner;
        val2 = dF->get_mat_knm(k, n, m);
        val3 = dLF->get_mat_knm(k, n, m) * ak;
        val4 = (dLH->get_mat_knm(k,n,m)+dLHN->get_mat_knm(k,n,m))*ak*besseli[n];
        set_mat_knm(k, n, m, val1+val2+val3+val4);
      }
    }
  }
}


GradFMatrix::GradFMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradFMatrix::calc_vals(Molecule mol, shared_ptr<IEMatrix> IE,
                            shared_ptr<GradWFMatrix> dWF)
{
  Ptx in, in2;
  for (int k = 0; k < get_ns(); k++)
  {
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
            in2=dWF->get_mat_knm(k,l,s)*IE->get_IE_k_nm_ls(k,n,m,l,s).real();
            in += in;
          }
        }
        set_mat_knm(k, n, m, in);
      }
    }
  }
}

GradHMatrix::GradHMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradHMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                            shared_ptr<IEMatrix> IE,
                            shared_ptr<GradWHMatrix> dWH)
{
  Ptx in, in2;
  double ak;
  vector<double> besseli;
  for (int k = 0; k < get_ns(); k++)
  {
    ak = mol.get_ak(k);
    besseli = bcalc->calc_mbfI(p_+1, ak * kappa_);
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
            in2 = dWH->get_mat_knm(k,l,s)*IE->get_IE_k_nm_ls(k,n,m,l,s).real();
            in += in2;
          }
        }
        set_mat_knm(k, n, m, in * besseli[n]);
      }
    }
  }
}

GradLFMatrix::GradLFMatrix(int I, int wrt, int ns, int p)
:GradCmplxMolMat(I, wrt, ns, p)
{
}

void GradLFMatrix::calc_vals(Molecule mol, shared_ptr<SHCalc> shcalc,
                             shared_ptr<TMatrix> T,
                             shared_ptr<GradFMatrix> dF, int Mp)
{
  MyMatrix<Ptx> reex, inner;
  for (int k = 0; k < get_ns(); k++)
  {
    inner = MyMatrix<Ptx> (p_, 2*p_+1);
    for (int j = 0; j < get_ns(); j++)
    {
      if (j == k) continue;
      if (T->is_analytic(I_, k, I_, j))
      {
        reex = T->re_expand_gradX(dF->get_mat_k(k), I_, k, I_, j);
      }
      else
      {
        reex = numeric_reex(mol, k, j, shcalc, dF);
      }
      inner += reex;
    }
    mat_[k] = inner;
  }
}

MyMatrix<Ptx> GradLFMatrix::numeric_reex(Molecule mol, int k, int j,
                                         shared_ptr<SHCalc> shcalc,
                                         shared_ptr<GradFMatrix> dF,
                                         int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, mol.get_ak(j));
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
        
        rb_k = sph_grid[b]-(mol.get_centerk(k)-mol.get_centerk(j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        
        inner = fpj * (1/rb_k.r());
        inner *= pow(mol.get_ak(k)/ rb_k.r(), n);
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

void GradLHMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                             shared_ptr<SHCalc> shcalc,
                             shared_ptr<TMatrix> T,
                             shared_ptr<GradHMatrix> dH, int Mp)
{
  MyMatrix<Ptx> reex, inner;
  
  vector<double> besseli, besselk;
  
  for (int k = 0; k < get_ns(); k++)
  {
    
    besseli = bcalc->calc_mbfI(p_+1, mol.get_ak(k) * kappa_);
    besselk = bcalc->calc_mbfK(p_+1, mol.get_ak(k) * kappa_);
    
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
        reex = numeric_reex(mol, k, j, besseli, besselk, shcalc, dH, Mp);
      }
      inner += reex;
    }
    mat_[k] = inner;
  }
}

MyMatrix<Ptx> GradLHMatrix::numeric_reex(Molecule mol, int k, int j,
                                         vector<double> besseli,
                                         vector<double> besselk,
                                         shared_ptr<SHCalc> shcalc,
                                         shared_ptr<GradHMatrix> dH,
                                         int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, mol.get_ak(j));
  Ptx hpj, inner;
  Pt rb_k;
  MyMatrix<Ptx> dLH_Ik (p_, 2 * p_ + 1);
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      for (int b = 0; b < Mp; b++)
      {
        hpj = calc_dh_P(sph_grid[b], k, besseli, shcalc, dH);
        hpj *= (4 * M_PI) / sph_grid.size();
        
        rb_k = sph_grid[b]-(mol.get_centerk(k)-mol.get_centerk(j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        
        inner = hpj * (1/rb_k.r());
        inner *= pow(mol.get_ak(k)/ rb_k.r(), n);
        inner *= conj(shcalc->get_result(n, m));
        inner *= exp(-kappa_*rb_k.r());
        inner *= besselk[n];
        
        dLH_Ik.set_val(n, m, inner);
      }
    }
  }
  return dLH_Ik;
}

Ptx GradLHMatrix::calc_dh_P(Pt P, int k, vector<double> besseli,
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
      in *= 1/besseli[n];
      dh += in;
    }
  }
  return dh;
}
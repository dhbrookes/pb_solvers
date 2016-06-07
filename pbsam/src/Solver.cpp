//
//  Solver.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solver.h"

ComplexMoleculeMatrix::ComplexMoleculeMatrix(int I, int ns, int p)
:p_(p), mat_(ns, MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
}

void ComplexMoleculeMatrix::reset_mat()
{
  for (int k = 0; k < mat_.size(); k++)
  {
    for (int i = 0; i < mat_[k].get_nrows(); i++)
    {
      for (int j = 0; j < mat_[k].get_ncols(); j++)
      {
        mat_[k].set_val(i, j, cmplx(0, 0));
      }
    }
  }
}

EMatrix::EMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void EMatrix::calc_vals(Molecule mol, shared_ptr<SHCalc> _shcalc, double eps_in)
{
  cmplx val;
  Pt cen;
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol.get_ns(); k++)
  {
    for (int alpha=0; alpha < mol.get_nc_k(k); alpha++)
    {
      cen = mol.get_posj(mol.get_ch_k_alpha(k, alpha));
      _shcalc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol.get_qj(mol.get_ch_k_alpha(k, alpha));
          r_alpha = mol.get_radj(mol.get_ch_k_alpha(k, alpha));
          a_k = mol.get_ak(k);
          
          val = conj(_shcalc->get_result(n, m));
          val *= q_alpha / eps_in;
          val *= pow(r_alpha / a_k, n);
          val += get_mat_knm(k, n, m);
          set_mat_knm(k, n, m, val);
        }
      }
    }
  }
}

LEMatrix::LEMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void LEMatrix::calc_vals(Molecule mol, shared_ptr<SHCalc> _shcalc,
                         double eps_in)
{
  cmplx val;
  Pt cen;
  
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol.get_ns(); k++)
  {
    for (int alpha=0; alpha < mol.get_nc_k(k); alpha++)
    {
      if (mol.get_cg_of_ch(alpha) == k) continue;
      cen = mol.get_posj_realspace(alpha) - mol.get_centerk(k);
      cen = mol.get_posj(mol.get_ch_k_alpha(k, alpha));
      _shcalc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol.get_qj(mol.get_ch_k_alpha(k, alpha));
          r_alpha = mol.get_radj(mol.get_ch_k_alpha(k, alpha));
          a_k = mol.get_ak(k);
          
          val = conj(_shcalc->get_result(n, m));
          val *= q_alpha / eps_in;
          val *= 1 / r_alpha;
          val *= pow(a_k / r_alpha, n);
          val += get_mat_knm(k, n, m);
          set_mat_knm(k, n, m, val);
        }
      }
    }
  }
}

IEMatrix::IEMatrix(int I, int ns, int p)
:p_(p), I_(I),
IE_(ns, MatOfMats<cmplx>::type(p, 2*p+1, MyMatrix<cmplx>(p, 2*p+1)))
{
}

IEMatrix::IEMatrix(int I, Molecule mol, shared_ptr<SHCalc> _shcalc, int p)
: p_(p), I_(I),
IE_(mol.get_ns(), MatOfMats<cmplx>::type(p, 2*p+1, MyMatrix<cmplx>(p, 2*p+1)))
{
  calc_vals(mol, _shcalc);
}

void IEMatrix::calc_vals(Molecule mol, shared_ptr<SHCalc> _shcalc)
{
  int m_grid = 20 * pow(p_, 2); //suggested in Yap 2010
  cmplx val;
  vector<Pt> grid_pts;
  Pt real_grid;
  for (int k = 0; k < mol.get_ns(); k++)
  {
    grid_pts = make_uniform_sph_grid(m_grid, mol.get_ak(k));
    // loop through other spheres in this molecule to find exposed grid points:
    for (int k2 = 0; k2< mol.get_ns(); k2++)
    {
      if (k2 == k) continue;
      for (int g = 0; g < m_grid; g++)
      {
        real_grid = grid_pts[g] - mol.get_centerk(k);
        // check if sphere is overlapping another
        if (real_grid.dist(mol.get_centerk(k2)) < mol.get_ak(k2)) continue;
        
        _shcalc->calc_sh(grid_pts[g].theta(), grid_pts[g].phi());
        for (int n = 0; n < p_; n++)
        {
          for (int m = -n; m < n+1; m++)
          {
            for (int l = 0; l < p_; l++)
            {
              for (int s = -l; s < l+1; s++)
              {
                val = get_IE_k_nm_ls(k, n, m, l, s);
                val += _shcalc->get_result(l, s) *
                conj(_shcalc->get_result(n, m));
                set_IE_k_nm_ls(k, n, m, l, s, val);
              }
            }
          }
        }
      }
    }
  }
}

TMatrix::TMatrix(int p, shared_ptr<System> _sys,
                 shared_ptr<SHCalc> _shcalc,
                 shared_ptr<Constants> _consts,
                 shared_ptr<BesselCalc> _besselcalc,
                 shared_ptr<ReExpCoeffsConstants> _reexpconsts)
:p_(p), kappa_(_consts->get_kappa()), Nmol_(_sys->get_n())
{
  int total_spheres=0;
  Nsi_ = vector<int> (Nmol_);
  for (int I = 0; I < Nmol_; I++)
  {
    total_spheres += _sys->get_Ns_i(I);
    Nsi_[I] = _sys->get_Ns_i(I);
  }
  T_.reserve(total_spheres*total_spheres);
  
  update_vals(_sys, _shcalc, _besselcalc, _reexpconsts);
}


void TMatrix::update_vals(shared_ptr<System> _sys, shared_ptr<SHCalc> _shcalc,
                          shared_ptr<BesselCalc> _besselcalc,
                          shared_ptr<ReExpCoeffsConstants> _reexpconsts)
{
  T_.clear();
  idxMap_.clear();
  
  int idx = 0;
  Pt c_Ik, c_Jl, v;
  vector<int> idx_vec;
  for (int I = 0; I < _sys->get_n(); I++)
  {
    for (int J = 0; J < _sys->get_n(); J++)
    {
      for (int k = 0; k < _sys->get_Ns_i(I); k++)
      {
        c_Ik = _sys->get_centerik(I, k);
        for (int l = 0; l < _sys->get_Ns_i(J); l++)
        {
          if (I==J && k==l) continue;
          
          idx_vec = {I, k, J, l};
          c_Jl = _sys->get_centerik(J, l);
          if (I==J && c_Ik.dist(c_Jl) < 5.0)
          {
            idxMap_[idx_vec] = -1;
            continue;
          }
          vector<double> besselK = _besselcalc->calc_mbfK(2*p_, kappa_ * v.r());
          v = _sys->get_pbc_dist_vec_base(c_Ik, c_Jl);
          _shcalc->calc_sh(v.theta(), v.phi());
          auto re_exp = make_shared<ReExpCoeffs>(p_, v,
                                                 _shcalc->get_full_result(),
                                                 besselK, _reexpconsts,
                                                 _sys->get_lambda(), false);
          T_.push_back(re_exp);
          idxMap_[idx_vec] = idx;
          idx++;
        }
      }
    }
  }
}

MyMatrix<cmplx> TMatrix::re_expand(int I, int k, int J, int l, MyMatrix<cmplx> X)
{
  MyMatrix<cmplx> X1 (p_, 2*p_ + 1);
  MyMatrix<cmplx> X2 (p_, 2*p_ + 1);
  MyMatrix<cmplx> Z (p_, 2*p_ + 1);
  
  int n, m, s;
  cmplx inter, sval, rval, aval;
  shared_ptr<ReExpCoeffs> T = get_T_Ik_Jl(I, k, J, l);
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (s = -n; s <= n; s++)
      {
          rval = T->get_rval(n, m, s);
          aval = X(n, s+p_);
          
          inter += rval * aval;
        } // end s
        X1.set_val(n, m+p_, inter);
      } // end m
    } //end n
  
  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        sval = T->get_sval(l, n, m);
        inter += sval * X1(l, m+p_);
      } // end l
      X2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  
  //fill z:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter = 0;
      for (s = -n; s <= n; s++)
      {
        rval = T->get_rval(n, s, m);
        inter += conj(rval) * X2(n, s+p_);
      } // end s
      Z.set_val(n, m+p_, inter);
    } // end m
  } //end n

  return Z;
}

LFMatrix::LFMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void LFMatrix::calc_vals(shared_ptr<TMatrix> T, shared_ptr<FMatrix> F,
                         shared_ptr<SHCalc> shcalc, shared_ptr<System> sys)
{
  reset_mat();
  MyMatrix<cmplx> reex;
  for (int k = 0; k < mat_.size(); k++)
  {
    for (int j = 0; j < T->get_nsi(I_); j++)
    {
      if (j==k) continue;
      if (! T->is_analytic(I_, k, I_, j))
      {
        reex = T->re_expand(I_, k, I_, j, F->get_mat_k(k));
      }
      else
      {
        reex = analytic_reex(I_, k, j, F, shcalc, sys);
      }
      mat_[k] += reex;
    }
  }
}

MyMatrix<cmplx> LFMatrix::analytic_reex(int I, int k, int j,
                                        shared_ptr<FMatrix> F,
                                        shared_ptr<SHCalc> shcalc,
                                        shared_ptr<System> sys, int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, sys->get_aik(I_, j));
  cmplx inner, fbj;
  Pt rb_k;
  MyMatrix<cmplx> LF_Ik (p_, 2 * p_ + 1);
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      for (int b = 0; b < Mp; b++)
      {
        fbj = make_fb_Ij(I, j, sph_grid[b], F, shcalc);
        fbj *= (4 * M_PI) / sph_grid.size();
        
        rb_k = sph_grid[b]-(sys->get_centerik(I, k)-sys->get_centerik(I, j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        
        inner = fbj / rb_k.r();
        inner *= pow(sys->get_aik(I, k) / rb_k.r(), n);
        inner *= conj(shcalc->get_result(n, m));
        LF_Ik.set_val(n, m, inner);
      }
    }
  }
  return LF_Ik;
}

cmplx LFMatrix::make_fb_Ij(int I, int j, Pt rb,
                                  shared_ptr<FMatrix> F,
                                  shared_ptr<SHCalc> shcalc)
{
  cmplx fbj, fb_inner;
  shcalc->calc_sh(rb.theta(),rb.phi());
  
  //create fbj function
  fbj = 0;
  for (int n2 = 0; n2 < p_; n2++)
  {
    for (int m2 = -n2; m2 < n2+1; m2++)
    {
      fb_inner = ((2 * n2 + 1) / (4 * M_PI)) * F->get_mat_knm(j, n2, m2);
      fb_inner *= shcalc->get_result(n2, m2);
      fbj += fb_inner;
    }
  }
  fbj /= pow(rb.r(), 2); // make f_hat into f
  return fbj;
}

LHMatrix::LHMatrix(int I, int ns, int p, double kappa)
:ComplexMoleculeMatrix(I, ns, p), kappa_(kappa)
{
}

void LHMatrix::calc_vals(shared_ptr<TMatrix> T, shared_ptr<HMatrix> H,
                         shared_ptr<SHCalc> shcalc, shared_ptr<System> sys,
                         shared_ptr<BesselCalc> bcalc, int Mp)
{
  reset_mat();
  MyMatrix<cmplx> reex;
  for (int k = 0; k < mat_.size(); k++)
  {
    for (int j = 0; j < T->get_nsi(I_); j++)
    {
      if (j==k) continue;
      if (! T->is_analytic(I_, k, I_, j))
      {
        reex = T->re_expand(I_, k, I_, j, H->get_mat_k(k));
      }
      else
      {
        reex = analytic_reex(I_, k, j, H, shcalc, sys, bcalc, Mp);
      }
      mat_[k] += reex;
    }
  }
}

MyMatrix<cmplx> LHMatrix::analytic_reex(int I, int k, int j,
                                        shared_ptr<HMatrix> H,
                                        shared_ptr<SHCalc> shcalc,
                                        shared_ptr<System> sys,
                                        shared_ptr<BesselCalc> bcalc,
                                        int Mp)
{
  if (Mp == -1) Mp = 2.5 * p_*p_;
  vector<Pt> sph_grid = make_uniform_sph_grid(Mp, sys->get_aik(I_, j));
  cmplx inner, hbj;
  Pt rb_k;
  MyMatrix<cmplx> LH_Ik (p_, 2 * p_ + 1);
  vector<double> besselk, besseli;

  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      for (int b = 0; b < Mp; b++)
      {
        
        hbj = make_hb_Ij(I, j, sph_grid[b], H, shcalc, bcalc);
        hbj *= (4 * M_PI) / (Mp);
        
        rb_k = sph_grid[b]-(sys->get_centerik(I, k)-sys->get_centerik(I, j));
        shcalc->calc_sh(rb_k.theta(), rb_k.phi());
        besselk = bcalc->calc_mbfK(p_+1, kappa_ * rb_k.r());
        
        inner = hbj / rb_k.r();
        inner *= pow(sys->get_aik(I, k) / rb_k.r(), n);
        inner *= exp(-kappa_*rb_k.r());
        inner *= besselk[n];
        inner *= conj(shcalc->get_result(n, m));
        LH_Ik.set_val(n, m, inner);
      }
    }
  }
  return LH_Ik;
}

cmplx LHMatrix::make_hb_Ij(int I, int j, Pt rb,
                           shared_ptr<HMatrix> H,
                           shared_ptr<SHCalc> shcalc,
                           shared_ptr<BesselCalc> bcalc)
{
  vector<double> besseli;
  cmplx hb_inner, hbj;
  vector<cmplx> h;

  shcalc->calc_sh(rb.theta(), rb.phi());
  besseli = bcalc->calc_mbfI(p_+1,
                             kappa_*rb.r());
  hbj = 0;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      hb_inner = ((2*n+1)/ (4*M_PI)) * H->get_mat_knm(j, n, m) ;
      hb_inner /= besseli[n];
      hb_inner *= shcalc->get_result(n, m);
      hbj += hb_inner;
    }
  }
  hbj /= pow(rb.r(), 2);  // make h_hat into h
  return hbj;
}

LHNMatrix::LHNMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void LHNMatrix::calc_vals(shared_ptr<TMatrix> T,
                          vector<shared_ptr<HMatrix> > H)
{
  reset_mat();
  MyMatrix<cmplx> reex;
  for (int k = 0; k < mat_.size(); k++)
  {
    for (int J = 0; J < T->get_nmol(); J++)
    {
      if (J == I_) continue;
      for (int l =0 ; l < T->get_nsi(J); l++)
      {
        reex = T->re_expand(I_, k, J, l, H[J]->get_mat_k(l));
        mat_[k] += reex;
      }
    }
  }
}

XHMatrix::XHMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void XHMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                         shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                         shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                         shared_ptr<LHNMatrix> LHN)
{
  cmplx inner;
  double ak;
  vector<double> in_k;
  for (int k = 0; k < mol.get_ns(); k++)
  {
    ak = mol.get_ak(k);
    in_k = bcalc->calc_mbfI(p_, k * ak);
    for (int n = 0; n < p_; n++)
    {
      for (int m = - n; m < n+1; m++)
      {
        inner = E->get_mat_knm(k, n, m);
        inner += ak*(LE->get_mat_knm(k, n, m)+LF->get_mat_knm(k, n, m));
        inner -= ak*in_k[n]*(LH->get_mat_knm(k, n, m)+LHN->get_mat_knm(k, n, m));
        set_mat_knm(k, n, m, inner);
      }
    }
  }
}


XFMatrix::XFMatrix(int I, int ns, int p, double eps_in, double eps_out)
:ComplexMoleculeMatrix(I, ns, p), eps_(eps_in/eps_out)
{
}

void XFMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                         shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                         shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                         shared_ptr<LHNMatrix> LHN, double kappa)
{
  cmplx inner;
  double ak;
  vector<double> in_k;
  for (int k = 0; k < mol.get_ns(); k++)
  {
    ak = mol.get_ak(k);
    in_k = bcalc->calc_mbfI(p_+1, k * ak);
    for (int n = 0; n < p_; n++)
    {
      for (int m = - n; m < n+1; m++)
      {
        inner = (pow(kappa * ak, 2) * in_k[n+1])/(2*n+3);
        inner += n * in_k[n];
        inner *= ak * (LH->get_mat_knm(k, n, m) + LHN->get_mat_knm(k, n, m));
        inner += (n+1) * eps_ * E->get_mat_knm(k, n, m);
        inner -= n*eps_*ak*(LE->get_mat_knm(k, n, m)+LF->get_mat_knm(k, n, m));
        set_mat_knm(k, n, m, inner);
      }
    }
  }
}

HMatrix::HMatrix(int I, int ns, int p, double kappa)
:ComplexMoleculeMatrix(I, ns, p), kappa_(kappa)
{
}


void HMatrix::calc_vals(Molecule mol, shared_ptr<HMatrix> prev,
                        shared_ptr<XHMatrix> XH,
                        shared_ptr<FMatrix> F,
                        shared_ptr<IEMatrix> IE,
                        shared_ptr<BesselCalc> bcalc)
{
  cmplx val, inner, inner2;
  double ak;
  vector<double> bessel_i, bessel_k;
  for (int k = 0; k < mol.get_ns(); k++)  // loop over each sphere
  {
    ak = mol.get_ak(k);
    bessel_i = bcalc->calc_mbfI(p_+1, kappa_*ak);
    bessel_k = bcalc->calc_mbfK(p_+1, kappa_*ak);
    for (int n = 0; n < p_; n++)  // rows in new matrix
    {
      for (int m = -n; m < n+1; m++)  // columns in new matrix
      {
        val = 0;
        for (int l = 0; l < p_; l++)  // rows in old matrix
        {
          for (int s = -l; s < l+1; s++)  //columns in old matrix
          {
            inner = (2 * l + 1) / bessel_i[l];
            inner -= exp(-kappa_*ak) * bessel_k[l];
            inner *= prev->get_mat_knm(k, l, s);
            inner += F->get_mat_knm(k, l, s);
            inner += XH->get_mat_knm(k, l, s);
            inner *= IE->get_IE_k_nm_ls(k, n, m, l, s);
            val += inner;
          }
        }
        inner *= bessel_i[n];
        set_mat_knm(k, n, m, val);
      }
    }
  }
}

FMatrix::FMatrix(int I, int ns, int p, double kappa)
:ComplexMoleculeMatrix(I, ns, p), kappa_(kappa)
{
}


void FMatrix::calc_vals(Molecule mol, shared_ptr<FMatrix> prev,
                        shared_ptr<XFMatrix> XF,
                        shared_ptr<HMatrix> H,
                        shared_ptr<IEMatrix> IE,
                        shared_ptr<BesselCalc> bcalc)
{
  cmplx val, inner, inner2;
  double ak;
  vector<double> bessel_i, bessel_k;
  for (int k = 0; k < mol.get_ns(); k++)  // loop over each sphere
  {
    ak = mol.get_ak(k);
    bessel_i = bcalc->calc_mbfI(p_+1, kappa_*ak);
    bessel_k = bcalc->calc_mbfK(p_+1, kappa_*ak);
    for (int n = 0; n < p_; n++)  // rows in new matrix
    {
      for (int m = -n; m < n+1; m++)  // columns in new matrix
      {
        val = 0;
        for (int l = 0; l < p_; l++)  // rows in old matrix
        {
          for (int s = -l; s < l+1; s++)  //columns in old matrix
          {
            inner = (l * bessel_k[l]) - ((2 * l +1) * bessel_k[l+1]);
            inner *= exp(-kappa_*ak) * H->get_mat_knm(k, l, s);
            inner += XF->get_mat_knm(k, l, s);
            inner += (2*l + 1 - l * XF->get_eps()) * prev->get_mat_knm(k, l, s);
            inner *= IE->get_IE_k_nm_ls(k, n, m, l, s);
            val += inner;
          }
        }
        set_mat_knm(k, n, m, val);
      }
    }
  }
}






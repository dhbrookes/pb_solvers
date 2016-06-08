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

// use Rakhmanov method
vector<Pt> IEMatrix::make_uniform_sph_grid(int m_grid, double r)
{
  
  vector<Pt> grid (m_grid);
  Pt gp;
  grid[0].set_r(r);
  grid[0].set_theta(0.0);
  grid[0].set_phi(0.0);
  double hk;
  for (int k = 0; k < m_grid; k++)
  {
    grid[k].set_r(r);
    hk = -1 + ((2 * (k-1)) / (m_grid-1));
    grid[k].set_theta(acos(hk));
    
    if (k==0 || k==m_grid-1) grid[k].set_phi(0);
    else
    {
      grid[k] = fmod(grid[k-1].phi() + (3.6/sqrt(m_grid) * (1/sqrt(1-(hk*hk)))),
                     2*M_PI);
    }
  }
  return grid;
}

TMatrix::TMatrix(int p, shared_ptr<System> _sys,
                 shared_ptr<SHCalc> _shcalc,
                 shared_ptr<Constants> _consts,
                 shared_ptr<BesselCalc> _besselcalc,
                 shared_ptr<ReExpCoeffsConstants> _reexpconsts)
:p_(p), kappa_(_consts->get_kappa())
{
  int total_spheres=0;
  for (int I = 0; I < _sys->get_n(); I++) total_spheres += _sys->get_Ns_i(I);
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
          if (c_Ik.dist(c_Jl) < 5.0)
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

MyMatrix<cmplx> TMatrix::re_expand(int I, int k, int j, int l, MyMatrix<cmplx> X)
{
  MyMatrix<cmplx> X1 (p_, 2*p_ + 1);
  MyMatrix<cmplx> X2 (p_, 2*p_ + 1);
  MyMatrix<cmplx> Z (p_, 2*p_ + 1);
  
  int n, m, s;
  cmplx inter, sval, rval, aval;
  shared_ptr<ReExpCoeffs> T = get_T_Ik_Jl(I, k, j, l);
  
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

LHMatrix::LHMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

LHNMatrix::LHNMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
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


XFMatrix::XFMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void XFMatrix::calc_vals(Molecule mol, shared_ptr<BesselCalc> bcalc,
                         shared_ptr<EMatrix> E, shared_ptr<LEMatrix> LE,
                         shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                         shared_ptr<LHNMatrix> LHN, double eps_in,
                         double eps_out, double kappa)
{
  cmplx inner;
  double ak;
  vector<double> in_k;
  double eps = eps_in / eps_out;
  for (int k = 0; k < mol.get_ns(); k++)
  {
    ak = mol.get_ak(k);
    in_k = bcalc->calc_mbfI(2*p_, k * ak);
    for (int n = 0; n < p_; n++)
    {
      for (int m = - n; m < n+1; m++)
      {
        inner = (pow(kappa * ak, 2) * in_k[n+1])/(2*n+3);
        inner += n * in_k[n];
        inner *= ak * (LH->get_mat_knm(k, n, m) + LHN->get_mat_knm(k, n, m));
        inner += (n+1) * eps * E->get_mat_knm(k, n, m);
        inner -= n*eps*ak*(LE->get_mat_knm(k, n, m)+LF->get_mat_knm(k, n, m));
        set_mat_knm(k, n, m, inner);
      }
    }
  }
}


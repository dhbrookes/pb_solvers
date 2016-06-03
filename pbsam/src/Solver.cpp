//
//  Solver.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solver.h"

EMatrix::EMatrix(int I, int ns, int p)
:p_(p), E_ (ns, MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
}

EMatrix::EMatrix(int i, Molecule mol, shared_ptr<SHCalc> _shcalc,
                 int p, double eps_in)
:p_(p), E_ (mol.get_ns(), MyMatrix<cmplx> (p, 2*p+1)), I_(i)
{
  calc_vals(mol, _shcalc, eps_in);
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
          val += get_E_knm(k, n, m);
          set_E_knm(k, n, m, val);
        }
      }
    }
  }
}

LEMatrix::LEMatrix(int I, int ns, int p)
:p_(p), LE_ (ns, MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
}

LEMatrix::LEMatrix(int I, Molecule mol, shared_ptr<SHCalc> _shcalc,
                 int p, double eps_in)
:p_(p), LE_ (mol.get_ns(), MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
  calc_vals(mol, _shcalc, eps_in);
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
          val += get_LE_knm(k, n, m);
          set_LE_knm(k, n, m, val);
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

vector<Pt> IEMatrix::make_uniform_sph_grid(int m_grid, double r)
{
  vector<Pt> grid (m_grid);
  Pt gp;
  double offset = 2.0/m_grid;
  double increment = M_PI * (3.0-sqrt(5.0));
  double x, y, z, rad, phi;
  for (int i = 0; i < m_grid; i++)
  {
    y = ((i * offset) - 1) + (offset / 2);
    rad = sqrt(1 - pow(y, 2));
    phi = (i % m_grid) * increment;
    x = cos(phi) * rad;
    z = sin(phi) * rad;

    grid[i].set_x(r*x);
    grid[i].set_y(r*y);
    grid[i].set_z(r*z);
  }
  return grid;
}


LFMatrix::LFMatrix(int I, int ns, int p)
:p_(p), LF_(ns, MyMatrix<cmplx>(p, 2*p+1)), I_(I)
{
}

LHMatrix::LHMatrix(int I, int ns, int p)
:p_(p), LH_(ns, MyMatrix<cmplx>(p, 2*p+1)), I_(I)
{
}

LHNMatrix::LHNMatrix(int I, int ns, int p)
:p_(p), LHN_(ns, MyMatrix<cmplx>(p, 2*p+1)), I_(I)
{
}


XHMatrix::XHMatrix(int I, int ns, int p)
:p_(p), XH_(ns, MyMatrix<cmplx>(p, 2*p+1)), I_(I)
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
        inner = E->get_E_knm(k, n, m);
        inner += ak*(LE->get_LE_knm(k, n, m)+LF->get_LF_knm(k, n, m));
        inner -= ak*in_k[n]*(LH->get_LH_knm(k, n, m)+LHN->get_LHN_knm(k, n, m));
        set_XH_knm(k, n, m, inner);
      }
    }
  }
}


XFMatrix::XFMatrix(int I, int ns, int p)
:p_(p), XF_(ns, MyMatrix<cmplx>(p, 2*p+1)), I_(I)
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
        inner *= ak * (LH->get_LH_knm(k, n, m) + LHN->get_LHN_knm(k, n, m));
        inner += (n+1) * eps * E->get_E_knm(k, n, m);
        inner -= n*eps*ak*(LE->get_LE_knm(k, n, m)+LF->get_LF_knm(k, n, m));
        set_XF_knm(k, n, m, inner);
      }
    }
  }
}


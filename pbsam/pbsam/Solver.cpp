//
//  Solver.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solver.h"

EMatrix::EMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc,
                 int p, double eps_in)
:p_(p), E_ (mol.get_ns(), MyMatrix<cmplx> (p, 2*p+1))
{
  cmplx val;
  Pt cen;
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol.get_ns(); k++)
  {
    for (int alpha=0; alpha < mol.get_nc_k(k); alpha++)
    {
      cen = mol.get_posj(mol.get_ch_k_alpha(k, alpha));
      sh_calc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol.get_qj(mol.get_ch_k_alpha(k, alpha));
          r_alpha = mol.get_radj(mol.get_ch_k_alpha(k, alpha));
          a_k = mol.get_ak(k);
          
          val = conj(sh_calc->get_result(n, m));
          val *= q_alpha / eps_in;
          val *= r_alpha / a_k;
          val += get_E_knm(k, n, m);
          set_E_knm(k, n, m, val);
        }
      }
    }
  }
}

LEMatrix::LEMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc,
                 int p, double eps_in)
:p_(p), LE_ (mol.get_ns(), MyMatrix<cmplx> (p, 2*p+1))
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
      sh_calc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol.get_qj(mol.get_ch_k_alpha(k, alpha));
          r_alpha = mol.get_radj(mol.get_ch_k_alpha(k, alpha));
          a_k = mol.get_ak(k);
          
          val = conj(sh_calc->get_result(n, m));
          val *= q_alpha / eps_in;
          val *= 1 / r_alpha;
          val *= a_k / r_alpha;
          val += get_LE_knm(k, n, m);
          set_LE_knm(k, n, m, val);
        }
      }
    }
  }
}


IEMatrix::IEMatrix(Molecule mol, shared_ptr<SHCalc> sh_calc, int p)
: p_(p),
IE_(mol.get_ns(), MatOfMats<cmplx>::type(p, 2*p+1, MyMatrix<cmplx>(p, 2*p+1)))
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
        
        sh_calc->calc_sh(grid_pts[g].theta(), grid_pts[g].phi());
        for (int n = 0; n < p_; n++)
        {
          for (int m = -n; m < n+1; m++)
          {
            for (int l = 0; l < p_; l++)
            {
              for (int s = -l; s < l+1; s++)
              {
                val = get_IE_k_nm_ls(k, n, m, l, s);
                val += sh_calc->get_result(l, s) *
                            conj(sh_calc->get_result(n, m));
                set_IE_k_nm_ls(k, n, m, l, s, val);
              }
            }
          }
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
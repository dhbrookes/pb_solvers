//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"

Molecule::Molecule(int M, double a, vector<double> qs, vector<Pt> pos, Pt cen)
:M_(M), a_(a), qs_(qs), pos_(pos), center_(cen)
{
  // repositioning the charges WRT center of charge
  for (int i = 0; i < M_; i++)
  {
    // check that the charge is encompassed by the the center and radius:
    if (pos_[i].dist(center_) > a_)
    {
      throw BadCenterException(center_, a_);
    }
    pos_[i] = pos_[i] - center_;
  }
}


Molecule::Molecule(int M, double a, vector<double> qs, vector<Pt> pos)
:M_(M), a_(a), qs_(qs), pos_(pos)
{
  // calculate the center of the molecule (for now using center):
  double xc, yc, zc;
  xc = 0;
  yc = 0;
  zc = 0;
  for (int i = 0; i < M_; i++)
  {
    xc += pos_[i].x();
    yc += pos_[i].y();
    zc += pos_[i].z();
  }
  xc /= M_;
  yc /= M_;
  zc /= M_;
  
  center_ = Pt(xc, yc, zc);
  
  // repositioning the charges WRT center of charge
  for (int i = 0; i < M_; i++)
  {
    pos_[i] = pos_[i] - center_;
  }
}


System::System(Constants consts, const vector<Molecule>& mols)
:consts_(consts), molecules_(mols), N_((int) mols.size())
{
  lambda_ = 0;
  for (int i = 0; i < N_; i++)
  {
    lambda_ += get_ai(i);
  }
  lambda_ /= N_;
}


const double System::calc_average_radius() const
{
  double ave = 0;
  for (int i = 0; i < N_; i++)
  {
    ave += get_ai(i);
  }
  ave  =  ave / N_;
  return ave;
}


void System::check_for_overlap()
{
  int i, j;
  double dist;
  Pt pi, pj;
  double ai, aj;
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      pi = molecules_[i].get_center();
      pj = molecules_[j].get_center();
      ai = molecules_[i].get_a();
      aj = molecules_[j].get_a();
      dist = pi.dist(pj);
      if (dist < (ai + aj)) throw OverlappingMoleculeException(i, j);
    }
  }
}

//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"


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
//    xc += abs(qs[i]) * pos_[i].x();
//    yc += abs(qs[i]) * pos_[i].y();
//    zc += abs(qs[i]) * pos_[i].z();
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
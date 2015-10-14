//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"


Molecule::Molecule(int M, double a, vector<double> qs, vector<EPt> pos)
:M_(M), a_(a), qs_(qs), pos_(pos)
{
}



System::System(Constants consts, const vector<Molecule>& mols)
:consts_(consts), molecules_(mols), N_((int) mols.size())
{
}


const double System::get_average_radius() const
{
  double ave = 0;
  for (int i = 0; i < N_; i++)
  {
    ave += get_ai(i);
  }
  ave  =  ave / N_;
  return ave;
}
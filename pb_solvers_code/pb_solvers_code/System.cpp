//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"

// user specified radius and center
Molecule::Molecule(string type, double a, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, double drot, double dtrans)
:type_(type), drot_(drot), dtrans_(dtrans), qs_(qs), pos_(pos), vdwr_(vdwr),
M_((int) pos.size()),
a_(a), center_(cen)
{
  reposition_charges();
}

// user specified radius
Molecule::Molecule(string type, double a, vector<double> qs, vector<Pt> pos,
         vector<double> vdwr, double drot, double dtrans)
:type_(type), drot_(drot), dtrans_(dtrans), qs_(qs), pos_(pos), vdwr_(vdwr),
M_((int) pos.size()),
a_(a)
{
  calc_center();
  reposition_charges();
}

// user specified center
Molecule::Molecule(string type, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, double drot, double dtrans)
:type_(type), drot_(drot), dtrans_(dtrans), qs_(qs), pos_(pos), vdwr_(vdwr),
M_((int) pos.size()),
center_(cen)
{
  reposition_charges();
  calc_a();
}

// neither the center or radius are specified
Molecule::Molecule(string type, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, double drot, double dtrans)
:type_(type), drot_(drot), dtrans_(dtrans), qs_(qs), pos_(pos), vdwr_(vdwr),
M_((int) pos.size())
{
  calc_center();
  reposition_charges();
  calc_a();
}


void Molecule::calc_center()
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
}

void Molecule::calc_a()
{
  a_ = 0;
  double dist;
  for (int i = 0; i < M_; i++)
  {
    dist = pos_[i].dist(center_) + vdwr_[i];
    if (dist > a_) a_ = dist;
  }
}

void Molecule::reposition_charges()
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


System::System(Constants consts, const vector<Molecule>& mols, double cutoff,
               double boxlength)
:consts_(consts), molecules_(mols), N_((int) mols.size()), cutoff_(cutoff),
boxLength_(boxlength)
{
  lambda_ = calc_average_radius();
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
    pi = molecules_[i].get_center();
    ai = molecules_[i].get_a();
    for (j = 0; j < N_; j++)
    {
      pj = molecules_[j].get_center();
      aj = molecules_[j].get_a();
      dist = pi.dist(pj);
      if (dist < (ai + aj)) throw OverlappingMoleculeException(i, j);
    }
  }
}

void Molecule::translate(Pt dr)
{
  center_ = center_ + dr;
}

void Molecule::rotate(Quat qrot)
{
  for (int i = 0; i < M_; i++)
  {
    cout << " This is new pt " << qrot.rotate_point(pos_[i]).x() << " " <<
    qrot.rotate_point(pos_[i]).y() << " " <<
    qrot.rotate_point(pos_[i]).z() << " " << endl;
    pos_[i] = qrot.rotate_point(pos_[i]);
  }
}


System::System(Constants consts, Setup setup, double cutoff)
:consts_(consts)
{
  vector<Molecule> mols;
  
  int i, j;
  string pqrpath;
  Molecule mol;
  for (i = 0; i < setup.get_ntype(); i++)
  {
    for (j = 0; j < setup.getTypeNCount(i); j++)
    {
      pqrpath = setup.getTypeNPQR(i)[j];
      PQRFile pqrobj (pqrpath);
      if (pqrobj.get_cg())  // coarse graining is in pqr
      {
        mol  = Molecule(setup.getTypeNDef(i), pqrobj.get_cg_radii()[0],
                        pqrobj.get_charges(),pqrobj.get_atom_pts(),
                        pqrobj.get_radii(), pqrobj.get_cg_centers()[0],
                        setup.getDrot(i), setup.getDtr(i));
      }
      else
      {
        mol = Molecule(setup.getTypeNDef(i), pqrobj.get_charges(),
                       pqrobj.get_atom_pts(), pqrobj.get_radii(),
                       setup.getDrot(i), setup.getDtr(i));
      }
      mols.push_back(mol);
    }
  }
  N_ = (int) mols.size();
  lambda_ = calc_average_radius();
  
  boxLength_ = setup.getBLen();
  cutoff_ = cutoff;
}

Pt System::get_pbc_dist_vec(int i, int j)
{
  Pt ci = get_centeri(i);
  Pt cj = get_centeri(j);
  Pt dv  = ci - cj;
  
  Pt v = Pt(dv.x() - round(dv.x()/boxLength_),
          dv.y() - round(dv.y()/boxLength_),
          dv.z() - round(dv.z()/boxLength_));

  return v;
}

bool System::less_than_cutoff(Pt v)
{
  if (v.norm() < cutoff_) return true;
  else return false;
}


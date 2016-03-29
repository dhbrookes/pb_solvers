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
}

// neither the center or radius are specified
Molecule::Molecule(string type, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, double drot, double dtrans)
:type_(type), drot_(drot), dtrans_(dtrans), qs_(qs), pos_(pos), vdwr_(vdwr),
M_((int) pos.size())
{
  calc_center();
  reposition_charges();
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
  xc /= (double) M_;
  yc /= (double) M_;
  zc /= (double) M_;
  
  center_ = Pt(xc, yc, zc);
}

void Molecule::calc_a()
{
  a_ = 0;
  double dist;
  for (int i = 0; i < M_; i++)
  {
    dist = pos_[i].norm() + vdwr_[i];
    if (dist > a_) a_ = dist;
  }
}

void Molecule::reposition_charges()
{
  bool recalc_a = false;
  // repositioning the charges WRT center of charge
  for (int i = 0; i < M_; i++)
  {
    // check that the charge is encompassed by the the center and radius:
    if (pos_[i].dist(center_) > a_)   recalc_a = true;
    pos_[i] = pos_[i] - center_;
  }
  
  if (recalc_a) calc_a();
}

void Molecule::translate(Pt dr)
{
  center_ = center_ + dr;
}

void Molecule::rotate(Quat qrot)
{
  for (int i = 0; i < M_; i++)
  {
    pos_[i] = qrot.rotate_point(pos_[i]);
  }
}

System::System(const vector<Molecule>& mols, double cutoff,
               double boxlength)
:molecules_(mols), N_((int) mols.size()), cutoff_(cutoff),
boxLength_(boxlength)
{
  check_for_overlap();
  lambda_ = calc_average_radius();
  if (boxLength_/2. < cutoff_)  compute_cutoff();
}

System::System(Setup setup, double cutoff)
:t_(0)
{
  vector<Molecule> mols;
  int i, j, chg;
  string pqrpath;
  Molecule mol;
  for (i = 0; i < setup.get_ntype(); i++)
  {
    PQRFile pqrI (setup.getTypeNPQR(i));
    XYZFile xyzI (setup.getTypeNXYZ(i), setup.getTypeNCount(i));
    for (j = 0; j < setup.getTypeNCount(i); j++)
    {
      vector<Pt> repos_charges(pqrI.get_M());
      Pt com = pqrI.get_cg_centers()[0];
      Pt move = xyzI.get_pts()[j] + com * -1.0;
      for ( chg = 0; chg < pqrI.get_M(); chg ++)
        repos_charges[chg] = pqrI.get_atom_pts()[chg] + move;
      if (pqrI.get_cg())  // coarse graining is in pqr
      {
        mol  = Molecule(setup.getTypeNDef(i), pqrI.get_cg_radii()[0],
                        pqrI.get_charges(), repos_charges,
                        pqrI.get_radii(), xyzI.get_pts()[j],
                        setup.getDrot(i), setup.getDtr(i));
      }
      else
      {
        mol = Molecule(setup.getTypeNDef(i), pqrI.get_charges(),
                       repos_charges, pqrI.get_radii(),
                       xyzI.get_pts()[j],
                       setup.getDrot(i), setup.getDtr(i));
      }
      molecules_.push_back(mol);
    }
  }
  N_ = (int) molecules_.size();
  boxLength_ = setup.getBLen();
  cutoff_ = cutoff;
  
  if (boxLength_/2. < cutoff_)  compute_cutoff();
  check_for_overlap();
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


void System::compute_cutoff()
{
  cutoff_ = boxLength_/2.0;
  cout << " The desired cutoff is larger than half the box length";
  cout << ". Resetting cutoff to 1/2 the boxlength: " << cutoff_ << endl;
}


void System::check_for_overlap()
{
  int i, j;
  double dist, ai, aj;
  for (i = 0; i < N_; i++)
  {
    ai = molecules_[i].get_a();
    for (j = i+1; j < N_; j++)
    {
      aj = molecules_[j].get_a();
      dist = get_pbc_dist_vec(i, j).norm();
      if (dist < (ai + aj)) throw OverlappingMoleculeException(i, j);
    }
  }
}

Pt System::get_pbc_dist_vec(int i, int j)
{
  Pt ci = get_centeri(i);
  Pt cj = get_centeri(j);
  Pt dv  = ci - cj;

  Pt v = Pt(dv.x() - round(dv.x()/boxLength_)*boxLength_,
          dv.y() - round(dv.y()/boxLength_)*boxLength_,
          dv.z() - round(dv.z()/boxLength_)*boxLength_);

  return v;
}

vector<Pt> System::get_allcenter() const
{
  vector< Pt> mol_cen(N_);
  for ( int i = 0; i < N_; i++)
    mol_cen[i] = molecules_[i].get_center();
  
  return mol_cen;
}

bool System::less_than_cutoff(Pt v)
{
  if (v.norm() < cutoff_) return true;
  else return false;
}

void System::write_to_pqr(string outfile)
{
  int i, j, ct = 0;
  ofstream pqr_out;
  char pqrlin[400];
  
  pqr_out.open( outfile );
  
  for ( i = 0; i < N_; i++ )
  {
    for ( j = 0; j < get_Mi(i); j++)
    {
      sprintf(pqrlin,"%6d  C   CHG A%-5d    %8.3f%8.3f%8.3f %7.4f %7.4f",ct,i,
              get_posijreal(i, j).x(),
              get_posijreal(i, j).y(),
              get_posijreal(i, j).z(),
              get_qij(i, j), get_radij(i, j));
      pqr_out << "ATOM " << pqrlin << endl;
      ct++;
    }
    sprintf(pqrlin,"%6d  X   CEN A%-5d    %8.3f%8.3f%8.3f %7.4f %7.4f",ct,i,
            get_centeri(i).x(), get_centeri(i).y(), get_centeri(i).z(),
            0.0, get_radi(i));
    pqr_out << "ATOM " << pqrlin << endl;
    ct++;
  }
  
}

//
//  System.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "System.h"


MoleculeCG::MoleculeCG(int type, int type_idx, string movetype, vector<double> qs,
          vector<Pt> pos, vector<double> vdwr, vector<Pt> cens,
          vector<double> as, double drot, double dtrans)
:type_(type), typeIdx_(type_idx), moveType_(movetype), qs_(qs), pos_(pos),
vdwr_(vdwr), centers_(cens), as_(as), drot_(drot), dtrans_(dtrans),
Nc_((int) qs.size()), Ns_((int) cens.size()), cgCharges_(cens.size())
{
  int closest;
  for (int i = 0; i < Nc_; i++)
  {
    closest = find_closest_center(pos_[i]);
    cgCharges_[closest].push_back(i);
    pos_[i] = pos_[i] - centers_[closest];  // reposition charge
    chToCG_[i] = closest;
  }
}

int MoleculeCG::find_closest_center(Pt pos)
{
  int idx = 0;
  double dmin = __DBL_MAX__;
  double d;
  for (int i = 0; i < Ns_; i++)
  {
    d = pos.dist(centers_[i]);
    if (d < dmin)
    {
      idx = i;
      dmin = d;
    }
  }
  return idx;
}

void MoleculeCG::translate(Pt dr, double boxlen)
{
  for (int k = 0; k < Ns_; k++)
  {
    Pt dv  = centers_[k] + dr;
    
    centers_[k] = Pt(dv.x() - round(dv.x()/boxlen)*boxlen,
                 dv.y() - round(dv.y()/boxlen)*boxlen,
                 dv.z() - round(dv.z()/boxlen)*boxlen);
  }
}

void MoleculeCG::rotate(Quat qrot)
{
  for (int k = 0; k < Ns_; k++)
  {
    centers_[k] = qrot.rotate_point(centers_[k]);
    for (int i = 0; i < cgCharges_[k].size(); i++)
    {
      pos_[cgCharges_[k][i]] = qrot.rotate_point(pos_[cgCharges_[k][i]]);
    }
  }
}

void MoleculeCG::find_centers(vector<Pt> sp, vector<Pt> np,
                              double tol_sp, int n_trials,
                              int max_trials)
{
  vector<int> unbound (Nc_);
  
  for (int i = 0; i < Nc_; i++) unbound[i] = i;
  int j = -1;
  int ct = 0;
  
  while(unbound.size() != 0 && ct < Nc_)
  {
    j++;
    centers_.resize(j+1);
    as_.resize(j+1);
    cgCharges_.resize(j+1);
    int n_max = 0;
    
    int m=-1;
    while (m < n_trials || (m >= n_trials && n_max==0 && m < max_trials))
    {
      m++;
      CGSphere best = find_best_center(sp, np, unbound, tol_sp);
      if (best.get_n() > n_max)
      {
        centers_[j] = best.get_center();
        as_[j] = best.get_a();
        cgCharges_[j] = best.get_ch();
        n_max = best.get_n();
      }
    }
    
    if (n_max == 0)
    {
      int rand_ind = unbound[(int) floor(drand48() * unbound.size())];
      cgCharges_[j] = { rand_ind };
      n_max = 1;
    }
    
    if (n_max == 1)
    {
      centers_[j] = pos_[cgCharges_[j][0]];
      as_[j] = vdwr_[cgCharges_[j][0]];
    }
    
    // remove bound charges from list of unbound charges:
    for (int l = 0; l < cgCharges_[j].size(); l++)
    {
      for (int f = 0; f < unbound.size(); f++)
      {
        if (unbound[f] == cgCharges_[j][l])
        {
          unbound.erase(unbound.begin() + f);
        }
      }
    }
    ct++;
  }
}

CGSphere MoleculeCG::find_best_center(vector<Pt> sp,vector<Pt> np,
                                      vector<int> unbound,
                                      double tol_sp, int iter, double beta)
{
  int sz = (int) unbound.size();
  const double sphere_tol =  1.0;
  
  Pt best_cen;
  double best_a = 0;
  int best_N = 0;
  vector<int> ch;  // encompassed charges of best sphere
  
  best_cen = pos_[unbound[(int) floor(drand48()*sz)]];
  
  for (int m = 0; m < iter; m++)
  {
    Pt tri_cen;
    double tri_a;
    int tri_N = 0;
    
    double cmin = __DBL_MAX__;
    
    for (int i = 0; i < sp.size(); i++)
    {
      double distsq = (sp[i] - best_cen).norm2();
      if (distsq < cmin) cmin = distsq;
    }
    double scale = sqrt(cmin);
    tri_cen = best_cen + (random_pt() * scale * random_normalized());
    
    int min_id = 0;
    tri_a = __DBL_MAX__;
    for (int i = 0; i < sp.size(); i++)
    {
      double distsq = (sp[i] - tri_cen).norm2();
      if (distsq < tri_a)
      {
        tri_a = distsq;
        min_id = i;
      }
    }
    if ((sp[min_id] - tri_cen).dot(np[min_id]) < 0.0) continue;
    tri_a = sqrt(tri_a) + tol_sp;
    tri_a *= tri_a;
    
    // count number of unbound charges within this sphere
    for (int i = 0; i < sz; i++)
    {
      double dist = (pos_[unbound[i]] - tri_cen).norm() + vdwr_[unbound[i]];
      if (dist < sqrt(tri_a)) tri_N++;
    }
    
    //apply MC criteria
    if (exp(beta*(tri_N - best_N) > drand48()))
    {
      best_cen = tri_cen;
      best_N = tri_N;
      best_a = tri_a;  // should this be sqrt() ??
    }
    
    if ((m+1) % 100 == 0) beta *= 1.1;
  }
  
  ch.reserve(best_N);
  for (int i = 0; i < sz; i++)
  {
    double dist = (pos_[unbound[i]] - best_cen).norm() + vdwr_[unbound[i]];
    if (dist < sqrt(best_a))
    {
      ch.push_back(unbound[i]);
    }
  }
  
  return CGSphere(best_cen, best_a, best_N, ch);
}




// user specified radius and center
Molecule::Molecule(string movetype, double a, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, int type, int typeIdx,
                   double drot, double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(typeIdx),
M_((int) pos.size()), a_(a), center_(cen)
{
  set_Dtr_Drot(movetype);
  reposition_charges();
}

// user specified radius
Molecule::Molecule(string movetype, double a, vector<double> qs,
                   vector<Pt> pos, vector<double> vdwr, int type, int typeIdx,
                   double drot, double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(typeIdx),
M_((int) pos.size()),
a_(a)
{
  set_Dtr_Drot(movetype);
  calc_center();
  reposition_charges();
}

// user specified center
Molecule::Molecule(string movetype, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, int type, int typeIdx,
                   double drot, double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(typeIdx),
M_((int) pos.size()), center_(cen), a_(0)
{
  set_Dtr_Drot(movetype);
  reposition_charges();
}

// neither the center or radius are specified
Molecule::Molecule(string movetype, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr,  int type, int typeIdx,
                   double drot, double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(typeIdx),
M_((int) pos.size()), a_(0)
{
  set_Dtr_Drot(movetype);
  calc_center();
  reposition_charges();
}


void Molecule::set_Dtr_Drot(string type)
{
  if ((type == "stat") or (type == "rot"))  dtrans_ = 0.0;
  if (type == "stat") drot_ = 0.0;
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
//  unwrappedCenter_ = center_;
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
    if (pos_[i].dist(center_)+vdwr_[i] > a_)   recalc_a = true;
    pos_[i] = pos_[i] - center_;
  }
  
  if (recalc_a) calc_a();
}

void Molecule::translate(Pt dr, double boxlen)
{
  Pt dv  = center_ + dr;
  
//  unwrappedCenter_ = unwrappedCenter_ + dr; // unwrapped position
  center_ = Pt(dv.x() - round(dv.x()/boxlen)*boxlen,
            dv.y() - round(dv.y()/boxlen)*boxlen,
            dv.z() - round(dv.z()/boxlen)*boxlen);
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
boxLength_(boxlength), t_(0)
{
  int i, j, k, maxi = 0;
  vector<int> maxj, keys(2);
  for ( k = 0; k < N_; k++)
  {
    i = molecules_[k].get_type();
    j = molecules_[k].get_type_idx();
    keys = {i,j};
    typeIdxToIdx_[keys] = k;
    maxi = ( maxi > i ) ? maxi : i;
    
    if ( i >= maxj.size() ) maxj.push_back(0);
    maxj[i] = ( maxj[i] > j ) ? maxj[i] : j;
  }
  
  maxi++;
  for ( j = 0; j < maxj.size(); j++) maxj[j]++;
  
  ntype_ = maxi;
  typect_ = maxj;
  
  check_for_overlap();
  lambda_ = calc_average_radius();
  if (boxLength_/2. < cutoff_)  compute_cutoff();
}

System::System(Setup setup, double cutoff)
:t_(0), ntype_(setup.get_ntype()), typect_(setup.get_type_nct())
{
  vector<Molecule> mols;
  int chg, i, j, k=0;
  string pqrpath;
  Molecule mol;
  vector<int> keys(2);
  for (i = 0; i < setup.get_ntype(); i++)
  {
    
    PQRFile pqrI (setup.getTypeNPQR(i));
    TransRotFile transrot;
    XYZFile xyzI;
    if (setup.getTypeIsTransRot(i))
    {
      transrot = TransRotFile(setup.getTypeNXYZ(i), setup.getTypeNCount(i));
    }
    else
    {
      xyzI = XYZFile (setup.getTypeNXYZ(i), setup.getTypeNCount(i));
    }
    for (j = 0; j < setup.getTypeNCount(i); j++)
    {
      Pt trans;
      MyMatrix<double> rot;
      
      Pt com = pqrI.get_cg_centers()[0];
      
      if (setup.getTypeIsTransRot(i))
      {
        trans = transrot.get_trans(j);
        rot = transrot.get_rotmat(j);
      }
      else
      {
        trans = xyzI.get_pts()[j] + com * -1.0;
        rot = MyMatrix<double> (3, 3, 0.0);
        rot.set_val(0, 0, 1.0);
        rot.set_val(1, 1, 1.0);
        rot.set_val(2, 2, 1.0);
      }
      
      keys = { i, j };
      vector<Pt> repos_charges(pqrI.get_M());
      Pt new_pt;
      
      for ( chg = 0; chg < pqrI.get_M(); chg ++)
      {
        repos_charges[chg] = pqrI.get_atom_pts()[chg].rotate(rot) + trans;
      }
      
      if (pqrI.get_cg())  // coarse graining is in pqr
      {
        mol  = Molecule(setup.getTypeNDef(i), pqrI.get_cg_radii()[0],
                        pqrI.get_charges(), repos_charges,
                        pqrI.get_radii(), xyzI.get_pts()[j], i, j,
                        setup.getDrot(i), setup.getDtr(i));
      }
      else if (! setup.getTypeIsTransRot(i))
      {
        mol = Molecule(setup.getTypeNDef(i), pqrI.get_charges(),
                       repos_charges, pqrI.get_radii(),
                       xyzI.get_pts()[j], i, j,
                       setup.getDrot(i), setup.getDtr(i));
        
      }
      else
      {
        mol = Molecule(setup.getTypeNDef(i), pqrI.get_charges(),
                       repos_charges, pqrI.get_radii(), i, j,
                       setup.getDrot(i), setup.getDtr(i));
      }
      
      molecules_.push_back(mol);
      typeIdxToIdx_[keys] = k;
      k++;
    } // end j
  } // end i
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
      if (dist < (ai + aj))
      {
        throw OverlappingMoleculeException(i, j);
      }
    }
  }
}

Pt System::get_pbc_dist_vec(int i, int j)
{
  Pt ci = get_centeri(i);
  Pt cj = get_centeri(j);
  return get_pbc_dist_vec_base(ci, cj);
}

Pt System::get_pbc_dist_vec_base(Pt p1, Pt p2)
{
  Pt dv  = p1 - p2;
  
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

void System::reset_positions( vector<string> xyzfiles )
{
  int i, j, k;
  vector<int> keys(2);
  for (i = 0; i < ntype_; i++)
  {
    XYZFile xyzI (xyzfiles[i], typect_[i]);
    for (j = 0; j < typect_[i]; j++)
    {
      keys = { i, j};
      k = typeIdxToIdx_[keys];
      Pt dist_to_new = get_centeri(k) - xyzI.get_pts()[j];
      molecules_[k].translate(dist_to_new*-1, boxLength_);
    }
  }
  
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

void System::write_to_xyz(ofstream & xyz_out)
{
  int i, j, at_tot = 0;
  char xyzlin[400];
  
  for ( i = 0; i < N_; i++ )
    for ( j = 0; j < get_Mi(i); j++)
      at_tot++;
  at_tot += N_; // for adding CG centers
  
  xyz_out << at_tot << endl;
  xyz_out << "Atoms. Timestep (ps): " << t_ << endl;
  for ( i = 0; i < N_; i++ )
  {
    for ( j = 0; j < get_Mi(i); j++)
    {
      sprintf(xyzlin,"N %8.3f %8.3f %8.3f", get_posijreal(i, j).x(),
              get_posijreal(i, j).y(), get_posijreal(i, j).z());
      xyz_out << xyzlin << endl;
    }
    sprintf(xyzlin,"X %8.3f %8.3f %8.3f", get_centeri(i).x(),
            get_centeri(i).y(), get_centeri(i).z());
    xyz_out << xyzlin << endl;
  }
}

//
//  SystemAM.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/28/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "SystemAM.h"

// user specified radius and center
MoleculeAM::MoleculeAM(string movetype, double a, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, int type, int typeIdx,
                   double drot, double dtrans)
:BaseMolecule(type, typeIdx, movetype, qs, pos, vdwr, cen, a, drot, dtrans),
unwrappedCenter_(cen)
{
  set_Dtr_Drot(movetype);
  reposition_charges();
}

// user specified radius
MoleculeAM::MoleculeAM(string movetype, double a, vector<double> qs,
                   vector<Pt> pos, vector<double> vdwr, int type, int typeIdx,
                   double drot, double dtrans)
:BaseMolecule(type, typeIdx, movetype, qs, pos, vdwr, drot, dtrans)
{
  set_Dtr_Drot(movetype);
  
  centers_[0] = calc_center();
  as_[0] = a;
  
  reposition_charges();
}

// user specified center
MoleculeAM::MoleculeAM(string movetype, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr, Pt cen, int type, int typeIdx,
                   double drot, double dtrans)
:BaseMolecule(type, typeIdx, movetype, qs, pos, vdwr, drot, dtrans)
{
  set_Dtr_Drot(movetype);
  centers_[0] = cen;
  as_[0] = 0;

  reposition_charges();
}

// neither the center or radius are specified
MoleculeAM::MoleculeAM(string movetype, vector<double> qs, vector<Pt> pos,
                   vector<double> vdwr,  int type, int typeIdx,
                   double drot, double dtrans)
:BaseMolecule(type, typeIdx, movetype, qs, pos, vdwr, drot, dtrans)
{
  set_Dtr_Drot(movetype);
  centers_[0] = calc_center();
  as_[0] = 0;
  reposition_charges();
}

Pt MoleculeAM::calc_center()
{
  // calculate the center of the MoleculeAM (for now using center):
  double xc(0), yc(0), zc(0);
  for (int i = 0; i < Nc_; i++)
  {
  //printf("Inside calc_center: cen: %.3f, %.3f, %.3f\n",
  //       pos_[i].x(), pos_[i].y(), pos_[i].z());
    xc += pos_[i].x();
    yc += pos_[i].y();
    zc += pos_[i].z();
  }
  xc /= (double) Nc_;
  yc /= (double) Nc_;
  zc /= (double) Nc_;
  
  Pt center = Pt(xc, yc, zc);
  return center;
}

double MoleculeAM::calc_a()
{
  double a = 0;
  double dist;
  for (int i = 0; i < Nc_; i++)
  {
    dist = pos_[i].norm() + vdwr_[i];
    if (dist > a) a = dist;
  }
  return a;
}

void MoleculeAM::reposition_charges()
{
  bool recalc_a = false;
  // repositioning the charges WRT center of charge
  for (int i = 0; i < Nc_; i++)
  {
    // check that the charge is encompassed by the the center and radius:
    if (pos_[i].dist(centers_[0])+vdwr_[i] > as_[0])
      recalc_a = true;
    pos_[i] = pos_[i] - centers_[0];
  //printf("Inside repos_charg: cen: %.3f, %.3f, %.3f\n",
  //       pos_[i].x(), pos_[i].y(), pos_[i].z());
  }
  
  if (recalc_a) as_[0] = calc_a();
}

void MoleculeAM::translate(Pt dr, double boxlen)
{
  Pt dv  = centers_[0] + dr;
  
  unwrappedCenter_ = unwrappedCenter_ + dr; // unwrapped position
  centers_[0] = Pt(dv.x() - round(dv.x()/boxlen)*boxlen,
            dv.y() - round(dv.y()/boxlen)*boxlen,
            dv.z() - round(dv.z()/boxlen)*boxlen);
}

void MoleculeAM::rotate(Quat qrot)
{
  for (int i = 0; i < Nc_; i++)
  {
    pos_[i] = qrot.rotate_point(pos_[i]);
  //cout << "Rotated pt " << pos_[i].x() << ", "
  //<< pos_[i].y() << ", "<< pos_[i].z() << endl;
  }
}


void MoleculeAM::rotate(MyMatrix<double> rotmat)
{
  for (int i = 0; i < Nc_; i++)
  {
    pos_[i] = pos_[i].rotate(rotmat);
  }
}

SystemAM::SystemAM(vector<shared_ptr<BaseMolecule> > mols, double cutoff,
               double boxlength)
:BaseSystem(mols, cutoff, boxlength)
{
  int i, j, k, maxi = 0;
  vector<int> maxj, keys(2);
  for ( k = 0; k < N_; k++)
  {
    i = molecules_[k]->get_type();
    j = molecules_[k]->get_type_idx();
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

SystemAM::SystemAM(Setup setup, double cutoff)
:BaseSystem()
{
  t_ = 0;
  ntype_ = setup.getNType();
  typect_ = setup.get_type_nct();
  
  vector<MoleculeAM> mols;
  int chg, i, j, k=0;
  string pqrpath;
  shared_ptr<MoleculeAM> mol;
  vector<int> keys(2);
  for (i = 0; i < setup.getNType(); i++)
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
      
      Pt com = pqrI.get_center_geo();
      
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
      vector<Pt> repos_charges(pqrI.get_Nc());
      Pt new_pt;
      
      for ( chg = 0; chg < pqrI.get_Nc(); chg ++)
      {
        repos_charges[chg] = pqrI.get_atom_pts()[chg].rotate(rot) + trans;
      }
      
      if (pqrI.get_Ns() != 0)  // coarse graining is in pqr
      {
        mol  = make_shared<MoleculeAM>(setup.getTypeNDef(i), pqrI.get_cg_radii()[0],
                        pqrI.get_charges(), repos_charges,
                        pqrI.get_radii(), xyzI.get_pts()[j], i, j,
                        setup.getDrot(i), setup.getDtr(i));
      }
      else if (! setup.getTypeIsTransRot(i))
      {
        mol = make_shared<MoleculeAM> (setup.getTypeNDef(i), pqrI.get_charges(),
                       repos_charges, pqrI.get_radii(),
                       xyzI.get_pts()[j], i, j,
                       setup.getDrot(i), setup.getDtr(i));
      }
      else
      {
        mol = make_shared<MoleculeAM > (setup.getTypeNDef(i), pqrI.get_charges(),
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

Pt SystemAM::get_pbc_dist_vec(int i, int j)
{
  Pt ci = get_centeri(i);
  Pt cj = get_centeri(j);
  return get_pbc_dist_vec_base(ci, cj);
}

vector<Pt> SystemAM::get_allcenter() const
{
  vector< Pt> mol_cen(N_);
  for ( int i = 0; i < N_; i++)
    mol_cen[i] = molecules_[i]->get_centerk(0);
  
  return mol_cen;
}

void SystemAM::reset_positions( vector<string> xyzfiles )
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
      molecules_[k]->translate(dist_to_new*-1, boxLength_);
    }
  }
  
}

void SystemAM::write_to_pqr(string outfile, int mid)
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

void SystemAM::write_to_xyz(ofstream & xyz_out)
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

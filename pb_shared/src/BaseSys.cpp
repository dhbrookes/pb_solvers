//
//  BaseSys.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 8/9/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BaseSys.h"



BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, double drot,
                           double dtrans)
:moveType_(movetype), drot_(drot), dtrans_(dtrans),
qs_(qs), pos_(pos), vdwr_(vdwr), type_(type), typeIdx_(type_idx),
Nc_((int) pos.size()), centers_(1), as_(1), Ns_(1)
{
  set_Dtr_Drot(moveType_);
}

BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, vector<Pt> cens,
                           vector<double> as, double drot, double dtrans)
:BaseMolecule(type, type_idx, movetype, qs, pos, vdwr, drot, dtrans)
{
  centers_ = cens;
  as_ = as;
  Ns_ = (int) centers_.size();
  
  set_Dtr_Drot(moveType_);
}

BaseMolecule::BaseMolecule(int type, int type_idx, string movetype,
                           vector<double> qs, vector<Pt> pos,
                           vector<double> vdwr, Pt cen,
                           double a, double drot, double dtrans)
:BaseMolecule(type, type_idx, movetype, qs, pos, vdwr, drot, dtrans)
{
  centers_[0] = cen;
  as_[0] = a;
  set_Dtr_Drot(moveType_);
}

void BaseMolecule::set_Dtr_Drot(string type)
{
  if ((type == "stat") or (type == "rot"))  dtrans_ = 0.0;
  if (type == "stat") drot_ = 0.0;
}


void BaseMolecule::write_pqr(string outfile)
{
  int j, k, ct(0);
  ofstream pqr_out;
  char pqrlin[400];
  
  pqr_out.open( outfile );
  
  for ( j = 0; j < get_nc(); j++)
  {
    sprintf(pqrlin,"%6d  C   CHG A%-5d    %8.3f%8.3f%8.3f %7.4f %7.4f",ct,0,
            get_posj_realspace(j).x(),
            get_posj_realspace(j).y(),
            get_posj_realspace(j).z(),
            get_qj(j), get_radj(j));
    pqr_out << "ATOM " << pqrlin << endl;
    ct++;
  }
  for (k = 0; k < get_ns(); k++)
  {
    sprintf(pqrlin,"%6d  X   CEN A%-5d    %8.3f%8.3f%8.3f %7.4f %7.4f",ct,0,
            get_centerk(k).x(), get_centerk(k).y(),
            get_centerk(k).z(), 0.0, get_ak(k));
    pqr_out << "ATOM " << pqrlin << endl;
    ct++;
  }
  pqr_out.close();
}

BaseSystem::BaseSystem(vector<shared_ptr<BaseMolecule> > mols,
                   double cutoff,
               double boxlength)
:molecules_(mols), N_((int) mols.size()), cutoff_(cutoff),
boxLength_(boxlength), t_(0)
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


const double BaseSystem::calc_average_radius() const
{
  double ave = 0;
  int total_sphere = 0;
  for (int i = 0; i < N_; i++)
  {
    total_sphere += get_Ns_i(i);
    for (int k = 0; k < get_Ns_i(i); k++)
    {
      ave += get_aik(i, k);
    }
  }
  ave  =  ave / (double) total_sphere;
  return ave;
}

void BaseSystem::compute_cutoff()
{
  cutoff_ = boxLength_/2.0;
  cout << " The desired cutoff is larger than half the box length";
  cout << ". Resetting cutoff to 1/2 the boxlength: " << cutoff_ << endl;
}

Pt BaseSystem::get_pbc_dist_vec_base(Pt p1, Pt p2)
{
  Pt dv  = p1 - p2;
  
  Pt v = Pt(dv.x() - round(dv.x()/boxLength_)*boxLength_,
            dv.y() - round(dv.y()/boxLength_)*boxLength_,
            dv.z() - round(dv.z()/boxLength_)*boxLength_);
  
  return v;
}

bool BaseSystem::less_than_cutoff(Pt v)
{
  if (v.norm() < cutoff_) return true;
  else return false;
}

void BaseSystem::check_for_overlap()
{
  int i, j, k1, k2;
  double dist, aik, ajk;
  Pt cen_ik, cen_jk;
  for (i = 0; i < N_; i++)
  {
    for (j = i+1; j < N_; j++)
    {
      for (k1 = 0; k1 < molecules_[i]->get_ns(); k1++)
      {
        for (k2 = 0; k2 < molecules_[j]->get_ns(); k2++)
        {
          aik = molecules_[i]->get_ak(k1);
          ajk = molecules_[j]->get_ak(k2);
          cen_ik = molecules_[i]->get_centerk(k1);
          cen_jk = molecules_[j]->get_centerk(k2);
          dist = get_pbc_dist_vec_base(cen_ik, cen_jk).norm();
          if (dist < (aik + ajk))
            throw OverlappingMoleculeException(i, j);
        }
      }
    }
  }
}



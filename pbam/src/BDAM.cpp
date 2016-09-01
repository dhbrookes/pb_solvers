//
//  BDAM.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BDAM.h"

ContactTerminateAM::ContactTerminateAM(vector<int> mol, double distance)
:BaseTerminate(), mol1_(mol[0]), mol2_(mol[1]), dist_contact_(distance)
{
  char buff[400];
  sprintf(buff, "Type %d and Type %d are within %5.2f;\t",
          mol1_, mol2_, dist_contact_);
  how_term_ = "System has fulfilled condition: " + string(buff);
}

const bool ContactTerminateAM::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  int i, j, idx1, idx2;
  for ( i = 0; i < _sys->get_typect(mol1_); i++)
  {
    for ( j = 0; j < _sys->get_typect(mol2_); j++)
    {
      idx1 = _sys->get_mol_global_idx( mol1_, i);
      idx2 = _sys->get_mol_global_idx( mol2_, j);
      
//      double mol_c2c = _sys->get_pbc_dist_vec(idx1, idx2).norm();
      double mol_c2c = _sys->get_pbc_dist_vec_base(_sys->get_centerik(idx1, 0),
                                                   _sys->get_centerik(idx2, 0)).norm();
      double a1 = _sys->get_aik(idx1, 0);
      double a2 = _sys->get_aik(idx2, 0);
      if (mol_c2c <= (dist_contact_+a1+a2)) return true;
    }
  }
  return false;
}


ContactTerminateAM2::ContactTerminateAM2(vector<int> mol,
                                         vector<vector<int> > atpairs,
                                         vector<double> dists, double pad)
:BaseTerminate(), mol1_(mol[0]), mol2_(mol[1]), atPairs_(atpairs),
dists_(dists), pad_(pad)
{
  string_create();
}

ContactTerminateAM2::ContactTerminateAM2(ContactFile confile, double pad)
:pad_(pad), mol1_(confile.get_moltype1()), mol2_(confile.get_moltype2()),
atPairs_(confile.get_at_pairs()), dists_(confile.get_dists())
{
  string_create();
}

void ContactTerminateAM2::string_create()
{
  char buff[400];
  sprintf(buff, "Type %d and Type %d are within %5.2f;\t",
          mol1_, mol2_, pad_);
  how_term_ = "System has fulfilled condition: " + string(buff);
}

const bool ContactTerminateAM2::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  bool contacted = false;
  int i, j, k, idx1, idx2;
  Pt cen1, cen2, pos1, pos2, vc1, vc2;
  double a1, a2, d, dcon;
  double sphdist1, sphdist2;  // distance of atom to edge of sphere
  
  for ( i = 0; i < _sys->get_typect(mol1_); i++)
  {
    for ( j = 0; j < _sys->get_typect(mol2_); j++)
    {
      idx1 = _sys->get_mol_global_idx( mol1_, i);
      idx2 = _sys->get_mol_global_idx( mol2_, j);
      
      cen1 = _sys->get_centerik(idx1, 0);
      cen2 = _sys->get_centerik(idx2, 0);
      
      a1 = _sys->get_aik(idx1, 0);
      a2 = _sys->get_aik(idx2, 0);
      
      for (k = 0; k < atPairs_.size(); k++)
      {
        dcon = dists_[k];
        pos1 = _sys->get_posijreal(idx1, atPairs_[k][0]);
        pos2 = _sys->get_posijreal(idx2, atPairs_[k][1]);
        
        vc1 = pos1 - cen1;
        vc2 = pos2 - cen2;
        
        sphdist1 = a1 - vc1.norm();
        sphdist2 = a2 - vc2.norm();
        
        // if sum of distances to edge of the spheres is > contact distance,
        // then contact can never happen and the new position is closest
        // point on edge of sphere and new contact distance is pad
        if ( (sphdist1 + sphdist2) > dcon)
        {
          pos1 = vc1 * (a1/vc1.norm()); // project onto sphere surface
          pos2 = vc2 * (a2/vc2.norm());
          dcon = pad_;
          
          // get position of atoms relative to box
          // (as opposed to center of MoleculeAM)
          pos1 = pos1 + cen1;
          pos2 = pos2 + cen2;
          d = _sys->get_pbc_dist_vec_base(pos1, pos2).norm();
          
          if (d < dcon){ contacted = true; break;}
        }else
        {
          d = _sys->get_pbc_dist_vec_base(pos1, pos2).norm();
          if (d < dcon)
          {
            contacted = true;
            cout << "This is dcon " << dcon << " and d " << d << endl;
            cout << "This is pos1 " << pos1.x() << " " <<
            pos1.y() << " " << pos1.z() << " " << endl;
            cout << "This is pos2 " << pos2.x() << " " <<
            pos2.y() << " " << pos2.z() << " " << endl;
            
            break;
          }
        }
      }
    }
  }
  return contacted;
}


BDStepAM::BDStepAM(shared_ptr<SystemAM> _sys, shared_ptr<Constants> _consts,
       vector<double> trans_diff_consts,
       vector<double> rot_diff_consts,
       bool diff, bool force)
:transDiffConsts_(trans_diff_consts), rotDiffConsts_(rot_diff_consts),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

BDStepAM::BDStepAM(shared_ptr<SystemAM> _sys, shared_ptr<Constants> _consts,
               bool diff, bool force)
:transDiffConsts_(_sys->get_n()), rotDiffConsts_(_sys->get_n()),
diff_(diff), force_(force), _sys_(_sys), _consts_(_consts)
{
  for (int i = 0; i < _sys_->get_n(); i++)
  {
    transDiffConsts_[i] = _sys_->get_dtransi(i);
    rotDiffConsts_[i] = _sys_->get_droti(i);
  }
  
  random_device rd;
  randGen_ = mt19937(rd());
}


double BDStepAM::compute_dt( )
{
  double DISTCUTOFF_TIME = 20.0;
  if ( min_dist_ - DISTCUTOFF_TIME > 0 )
    return 2.0 + ( min_dist_ - DISTCUTOFF_TIME )/15.0;
  else
    return 2.0;
}

void BDStepAM::compute_min_dist( )
{
  int i, j;
  Pt pt1, pt2;
  double newD, minDist = 1e100;
  for (i = 0; i < _sys_->get_n(); i++)
  {
    for (j = i+1; j < _sys_->get_n(); j++)
    {
      newD = _sys_->get_pbc_dist_vec(i, j).norm();
      if (newD < minDist) minDist = newD;
    }
  }
  
  min_dist_ = minDist;
}

Pt BDStepAM::rand_vec(double mean, double var)
{
  normal_distribution<double> dist (mean, sqrt(var));
  Pt pout = Pt(dist(randGen_), dist(randGen_), dist(randGen_));
  return pout;
}

void BDStepAM::indi_trans_update(int i, Pt fi)
{
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = transDiffConsts_[i] * dt_ * ikT_int;
  Pt dr = Pt(fi * coeff);
  Pt rand, new_pt;
  bool accept = false;
  
  int tries = 0;
  while (!accept && tries < 500)
  {
    tries++;
    rand = (diff_) ? rand_vec(0,2*transDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    _sys_->translate_mol(i, dr + rand);
    accept = true;
    try {
      _sys_->check_for_overlap();
    } catch (OverlappingMoleculeException) {
      _sys_->translate_mol(i, (dr+rand) * -1);
      accept = false;
    }
  }
}


void BDStepAM::indi_rot_update(int i, Pt tau_i)
{
  Pt rand, new_pt;
  Quat qrot;
  
  bool accept = false;
  double kT = _consts_->get_kbt();
  double ikT_int = 1 / Constants::convert_j_to_int(kT);
  double coeff = rotDiffConsts_[i] * dt_ * ikT_int;
  
  Pt dtheta = (tau_i * coeff);
  Pt center = _sys_->get_centeri(i);
  
  while (! accept)
  {
    rand = (diff_) ? rand_vec(0,2*rotDiffConsts_[i]*dt_) : Pt(0.0,0.0,0.0);
    dtheta = dtheta + rand;
    
    // creating zero quaternion if there is no rot
    if (dtheta.norm() < 1e-15) qrot = Quat();
    else qrot = Quat(dtheta.norm(), dtheta);
//    new_pt = qrot.rotate_point(center);
    _sys_->rotate_mol(i, qrot);
    accept = true;
    try {
      _sys_->check_for_overlap();
    } catch (OverlappingMoleculeException) {
      _sys_->rotate_mol(i, qrot.conj());
      accept = false;
    }
//    if (!check_for_collision(i, new_pt))
//    {
//      _sys_->rotate_mol(i, qrot);
//      accept = true;
//    }
  }
}

void BDStepAM::bd_update(shared_ptr<vector<Pt> > _F,
                   shared_ptr<vector<Pt> > _tau)
{
  int i;
  compute_min_dist();
  dt_ = compute_dt();
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( transDiffConsts_[i] != 0) indi_trans_update(i, _F->operator[](i));
  }
  
  for (i = 0; i < _sys_->get_n(); i++)
  {
    if ( rotDiffConsts_[i] != 0) indi_rot_update(i, _tau->operator[](i));
  }
  update_sys_time(dt_);
}


BDRunAM::BDRunAM(shared_ptr<ASolver> _asolv,
             shared_ptr<BaseTerminate> _terminator, string outfname, int num,
             bool diff, bool force, int maxiter, double prec)
:maxIter_(maxiter), _asolver_(_asolv), prec_(prec), _terminator_(_terminator)
{
  if (num == 0) _physCalc_ = make_shared<PhysCalcAM>(_asolver_, outfname);
  else _physCalc_ = make_shared<ThreeBodyPhysCalcAM>(_asolver_, num, outfname);
  
  _stepper_ = make_shared<BDStepAM> (_asolver_->get_sys(),
                                   _asolver_->get_consts(), diff, force);
}

void BDRunAM::run(string xyzfile, string statfile, int nSCF)
{
  int i(0), scf(2);
  int WRITEFREQ = 200;
  bool term(false), polz(true);
  ofstream xyz_out, stats;
  xyz_out.open(xyzfile);
  stats.open(statfile, fstream::in | fstream::out | fstream::app);
  
  while (i < maxIter_ and !term)
  {
    if ((i % WRITEFREQ) == 0 )
    {
      _stepper_->get_system()->write_to_xyz(xyz_out);
      if (i != 0)  _physCalc_->print_all();
    }
    
    _stepper_->get_system()->clear_all_lists();
    _asolver_->reset_all(_stepper_->get_system());
    if (nSCF != 0) scf = nSCF;

    _asolver_->solve_A(prec_, scf);
    _asolver_->solve_gradA(prec_, scf);
    _physCalc_->calc_force();
    _physCalc_->calc_torque();
    _stepper_->bd_update(_physCalc_->get_F(), _physCalc_->get_Tau());
   
    /*int N = _stepper_->get_system()->get_n();
      auto fo = make_shared<vector<Pt> > (N, Pt(0.0,0.0,0.0));
      auto to = make_shared<vector<Pt> > (N, Pt(0.0,0.0,0.0));
      _stepper_->bd_update(fo,to);
   */

    if ( (i % 100) == 0 ) cout << "This is step " << i << " and polz " << polz<< endl;

    if (_terminator_->is_terminated(_stepper_->get_system()))
    {
      term = true;
      // Printing out details at end
      _stepper_->get_system()->write_to_xyz(xyz_out);
      _physCalc_->print_all();
      stats << _terminator_->get_how_term(_stepper_->get_system());
      stats << " at time (ps) " << _stepper_->get_system()->get_time() << endl;
    }
    i++;
  }
  
  if ( i >= maxIter_ )
    stats << "System has gone over max number of BD iterations" << endl;
  
  xyz_out.close();
  stats.close();
}

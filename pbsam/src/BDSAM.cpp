//
//  BDSAM.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 1/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "BDSAM.h"

ContactTerminateSAM::ContactTerminateSAM(vector<int> mol,
                                         vector<vector<int> > atpairs,
                                         vector<double> dists, double pad)
:BaseTerminate(), mol1_(mol[0]), mol2_(mol[1]), atPairs_(atpairs),
dists_(dists), pad_(pad)
{
  string_create();
}


ContactTerminateSAM::ContactTerminateSAM(ContactFile confile, double pad)
:pad_(pad), mol1_(confile.get_moltype1()), mol2_(confile.get_moltype2()),
atPairs_(confile.get_at_pairs()), dists_(confile.get_dists())
{
  string_create();
}

void ContactTerminateSAM::string_create()
{
  char buff[400];
  sprintf(buff, "Type %d and Type %d are in contact;\t",
          mol1_, mol2_);
  how_term_ = "System has fulfilled condition: " + string(buff);
}

const bool ContactTerminateSAM::is_terminated(shared_ptr<BaseSystem> _sys) const
{
  bool contacted = false;
  int i, j, ctct, k, idx1, idx2, sph1, sph2;
  Pt cen1, cen2, vc1, vc2;
  double a1, a2, dcon;
  
  for ( i = 0; i < _sys->get_typect(mol1_); i++)
  {
    for ( j = 0; j < _sys->get_typect(mol2_); j++)
    {
      ctct = 0;
      idx1 = _sys->get_mol_global_idx( mol1_, i);
      idx2 = _sys->get_mol_global_idx( mol2_, j);
      
      for (k = 0; k < atPairs_.size(); k++)
      {
        dcon = dists_[k];
        
        sph1 = _sys->get_moli(idx1)->get_cg_of_ch(atPairs_[k][0]);
        sph2 = _sys->get_moli(idx1)->get_cg_of_ch(atPairs_[k][0]);
        
        cen1 = _sys->get_centerik(idx1, sph1);
        cen2 = _sys->get_centerik(idx2, sph2);
        
        a1 = _sys->get_aik(idx1, sph1);
        a2 = _sys->get_aik(idx2, sph2);
        
        // if the distance between the 2 CG spheres is less than the cutoff
        // The points are in contact
        if ( _sys->get_pbc_dist_vec_base( cen1, cen2).norm() < a1+a2+dcon)
        {
          cout << "This is dcon " << dcon << " sph1 " << cen1.x() << " " <<
          cen1.y() << " " << cen1.z() << " " << " & sph2 " << cen2.x()
          << " " << cen2.y() << " " << cen2.z() << " " << endl;
          ctct++;
          if (ctct == 2){contacted = true; break;}
        }
      }
    }
  }
  return contacted;
}

BDStepSAM::BDStepSAM(shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
       vector<double> trans_diff_consts,
       vector<double> rot_diff_consts,
       bool diff, bool force)
:BaseBDStep(_sys, _consts, trans_diff_consts, rot_diff_consts, diff, force)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

BDStepSAM::BDStepSAM(shared_ptr<BaseSystem> _sys, shared_ptr<Constants> _consts,
               bool diff, bool force)
:BaseBDStep(_sys, _consts, diff, force)
{
  random_device rd;
  randGen_ = mt19937(rd());
}

void BDStepSAM::compute_min_dist( )
{
  
  int i, j;
  Pt pt1, pt2;
  double newD, minDist = 1e100;
  for (i = 0; i < _sys_->get_n(); i++)
  {
    for (j = i+1; j < _sys_->get_n(); j++)
    {
      newD = _sys_->calc_min_dist(i, j);
      if (newD < minDist) minDist = newD;
    }
  }
  min_dist_ = minDist;
}



BDRunSAM::BDRunSAM(shared_ptr<Solver> _solv, shared_ptr<GradSolver> _gradSolv,
             shared_ptr<BaseTerminate> _terminator, string outfname,
             bool diff, bool force, int maxiter, double prec)
:BaseBDRun(_terminator, outfname, 0, diff, force, maxiter, prec),
_solver_(_solv), _gradSolv_(_gradSolv)
{
  _physCalc_ = make_shared<PhysCalcSAM>(_solv, _gradSolv, outfname);
  _stepper_ = make_shared<BDStepSAM> (_solver_->get_sys(),
                                   _solver_->get_consts(), diff, force);
}

void BDRunSAM::run(string xyzfile, string statfile, int nSCF)
{
  int i(0), scf(2), WRITEFREQ(200);
  bool term(false);
  ofstream xyz_out, stats;
  xyz_out.open(xyzfile);
  stats.open(statfile, fstream::in | fstream::out | fstream::app);

  while (i < maxIter_ and !term)
  {
    _solver_->reset_all();
    if (nSCF != 0) scf = nSCF;

    _solver_->solve(prec_, scf);
    _gradSolv_->update_HF(_solver_->get_all_F(), _solver_->get_all_H());
    _gradSolv_->solve(prec_, scf);
    _physCalc_->calc_force();
    _physCalc_->calc_torque();
    
    
    if ((i % WRITEFREQ) == 0 )
    {
      _stepper_->get_system()->write_to_xyz(xyz_out);
//      cout << "This is step " << i << endl;
//      for (int i = 0; i<_stepper_->get_system()->get_n(); i++)
//      {
//      cout << "This is force " <<  _physCalc_->get_forcei(i).x() <<", "<< _physCalc_->get_forcei(i).y() << ", " << _physCalc_->get_forcei(i).z() <<  endl;
//      }
        _physCalc_->print_all();
    }
    
    _stepper_->bd_update(_physCalc_->get_F(), _physCalc_->get_Tau());
    
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

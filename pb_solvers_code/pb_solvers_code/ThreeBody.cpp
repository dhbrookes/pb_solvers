//
//  ThreeBody.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/9/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "ThreeBody.h"

ThreeBody::ThreeBody( shared_ptr<ASolver> _asolver, Units unt, double cutoff )
: N_(_asolver->get_N()), p_(_asolver->get_p()), cutoffTBD_(cutoff),
_besselCalc_(_asolver->get_bessel()),
_shCalc_(_asolver->get_sh()),
_consts_(_asolver->get_consts()),
_sys_(_asolver->get_sys()),
energy_approx_(N_, 0),
force_approx_(N_, Pt(0,0,0)),
torque_approx_(N_, Pt(0,0,0)),
unit_(unt)
{
  dimer_.reserve(N_*N_);
  trimer_.reserve(N_*N_*N_);
  
  generatePairsTrips();
}

void ThreeBody::generatePairsTrips()
{
  int i, j, k;
  vector<double> dist(3); // distances between pairs: [ij, ik, jk]
  vector<int> temp2(2), temp3(3);
  int M2 = N_*(N_-1)/2; int M3 = N_*(N_-1)*(N_-2)/6;
  
  for( i = 0; i < N_; i++)
  {
    for( j = i+1; j < N_; j++)
    {
      dist[0] = _sys_->get_pbc_dist_vec( i, j).norm();
      if (cutoffTBD_ > dist[0])
      {
        temp2[0] = i; temp2[1] = j;
        dimer_.push_back(temp2);
      }
      for( k = j + 1; k < N_; k++)
      {
        dist[1] = _sys_->get_pbc_dist_vec( i, k).norm();
        dist[2] = _sys_->get_pbc_dist_vec( j, k).norm();
        sort( dist.begin(), dist.end());
        if (cutoffTBD_*2.0 > (dist[0] + dist[1]))
        {
          temp3[0] = i; temp3[1] = j; temp3[2] = k;
          trimer_.push_back(temp3);
        }
      }
    }
  }
  
  // Resizing vectors for saving energy/force vals
  energy_di_.resize( dimer_.size());
  force_di_.resize( dimer_.size());
  torque_di_.resize( dimer_.size());
  for ( i = 0; i < dimer_.size(); i++)
  {
    energy_di_[i].resize(2);
    force_di_[i].resize(2);
    torque_di_[i].resize(2);
  }
  
  energy_tri_.resize( trimer_.size());
  force_tri_.resize( trimer_.size());
  torque_tri_.resize( trimer_.size());
  for ( i = 0; i < trimer_.size(); i++)
  {
    energy_tri_[i].resize(3);
    force_tri_[i].resize(3);
    torque_tri_[i].resize(3);
  }
  
  if ( cutoffTBD_ < 1e37 )
    cout << "Cutoffs implemented, using " << cutoffTBD_ << endl;
  cout << "Max # of di: "  << M2 << " Act used: " << dimer_.size() <<endl;
  cout << "Max # of tri: " << M3 << " Act used: " << trimer_.size() <<endl;
} //end cutoffTBD


shared_ptr<System> ThreeBody::make_subsystem(vector<int> mol_idx)
{
  vector<Molecule> sub_mols (mol_idx.size());
  for (int i = 0; i < mol_idx.size(); i++)
  {
    sub_mols[i] = _sys_->get_molecule(mol_idx[i]);
  }
  
  shared_ptr<System> _subsys = make_shared<System>(sub_mols,_sys_->get_cutoff(),
                                                  _sys_->get_boxlength());
  _subsys -> set_time(_sys_->get_time());
  return _subsys;
}

// Two or three body approximation computation
void ThreeBody::solveNmer( int num, double preclim )
{
  int i, j;
  shared_ptr<vector<vector<int> > > nmer = (( num == 2 ) ?
                              make_shared<vector<vector<int> > >(dimer_) :
                              make_shared<vector<vector<int> > >(trimer_));
  vector< Molecule > mol_temp;
  shared_ptr<System> _sysTemp = make_subsystem(nmer->operator[](0));
  shared_ptr<ASolver> _asolvTemp = make_shared<ASolver>(_besselCalc_, _shCalc_,
                                                        _sysTemp, _consts_, p_);
  
  for( i = 0; i < nmer->size(); i++)
  {
    shared_ptr<System> _sysTemp = make_subsystem(nmer->operator[](i));
    
    _asolvTemp->reset_all(_sysTemp);
    _asolvTemp->solve_A(preclim);
    _asolvTemp->solve_gradA(preclim);
    
    PhysCalc phys_all( _asolvTemp, unit_);
    phys_all.calc_all();
    
    for ( j = 0; j < num; j++)
    {
      if ( num == 2 )
      {
        energy_di_[i][j] = phys_all.get_omegai_conv(j);
        force_di_[i][j] = phys_all.get_forcei_conv(j);
        torque_di_[i][j] = phys_all.get_taui_conv(j);
      } else
      {
        energy_tri_[i][j] = phys_all.get_omegai_conv(j);
        force_tri_[i][j] = phys_all.get_forcei_conv(j);
        torque_tri_[i][j] = phys_all.get_taui_conv(j);
      }
    }
  }
  
  cout << num << "mers done " << endl;
}

// Three body approximation
int ThreeBody::find_di( int i, int j)
{
  int di;
  vector<int> dim(2);
  
  for( di = 0; di < dimer_.size(); di++)
  {
    if (( dimer_[di][0] == i ) and ( dimer_[di][1] == j ))
      return di;
  }
  
  return -1;
}

// Three body approximation
void ThreeBody::calcTBDEnForTor( )
{
  int i, j, k;
  vector<int> m(3), di(3); // Given any triplet, di = [01, 02, 12] pairs
  // Matrix of if 0, 1 or 2 mol is in the [01, 02, 12] pairs & loc if so
  vector<vector<int > > loc(3, vector<int> (3));
  vector<vector<double > > incl(3, vector<double> (3));
  incl[0][0] = 1; incl[0][1] = 1; incl[0][2] = 0;
  incl[1][0] = 1; incl[1][1] = 0; incl[1][2] = 1;
  incl[2][0] = 0; incl[2][1] = 1; incl[2][2] = 1;
  
  loc[0][0] = 0; loc[0][1] = 0; loc[0][2] = 0;
  loc[1][0] = 1; loc[1][1] = 0; loc[1][2] = 0;
  loc[2][0] = 0; loc[2][1] = 1; loc[2][2] = 1;
  
  calcTwoBDEnForTor( ); // 2BD contribution
  
  for( i = 0; i < trimer_.size(); i++)
  {
    for( j = 0; j < 3; j++) m[j] = trimer_[i][j];
    
    di[0] = find_di(m[0], m[1]);
    di[1] = find_di(m[0], m[2]);
    di[2] = find_di(m[1], m[2]);
    
    for( j = 0; j < 3; j++)
    {
      energy_approx_[m[j]] += energy_tri_[i][j];
      force_approx_[m[j]]   = force_approx_[m[j]] + force_tri_[i][j];
      torque_approx_[m[j]]  = torque_approx_[m[j]] + torque_tri_[i][j];
      
      for ( k = 0; k < 3; k++)
      {
        energy_approx_[m[j]] -= incl[j][k]*energy_di_[di[k]][loc[j][k]];
        force_approx_[m[j]] =  (force_approx_[m[j]] +
                                force_di_[di[k]][loc[j][k]]*incl[j][k]*-1);
        torque_approx_[m[j]] =  (torque_approx_[m[j]] +
                                 torque_di_[di[k]][loc[j][k]]*incl[j][k]*-1);
        
      } // end k
    } // end j
  } // end i
}

void ThreeBody::printTBDEnForTor( vector<string> outfile )
{
  int j;
  for ( j = 0; j < N_; j++)
  {
    cout << "This is mol " << j << endl;
    cout << " Energy: " << get_energyi_approx(j) << "\t";
    cout << " Force: " << get_forcei_approx(j).norm() << "\t [";
    cout << get_forcei_approx(j).x();
    cout << ", " << get_forcei_approx(j).y()<< ", ";
    cout << get_forcei_approx(j).z()<< "]"<< endl;
  }
  
  if ( outfile[0] != "") printNmer( 2, outfile[0]);
  if ( outfile[1] != "") printNmer( 3, outfile[1]);
}

void ThreeBody::printNmer( int num, string outfile)
{
  int i, j;
  ofstream nmer_deets;
  double dist, en_nrm;
  Pt fo_nrm;
  vector< double > print_all(4);
  vector<int>  mol(3);
  shared_ptr<vector<vector<int> > > nmer = (( num == 2 ) ?
                                            make_shared<vector<vector<int> > >(dimer_) :
                                            make_shared<vector<vector<int> > >(trimer_));
  
  shared_ptr<vector<vector<double > > > en = ( num == 2 ) ?
  make_shared<vector<vector<double> > >(energy_di_) :
  make_shared<vector<vector<double > > >(energy_tri_);
  shared_ptr<vector<vector<Pt> > > frc = (( num == 2 ) ?
                                          make_shared<vector<vector<Pt> > >(force_di_) :
                                          make_shared<vector<vector<Pt> > >(force_tri_));
  shared_ptr<vector<vector<Pt> > > tor = (( num == 2 ) ?
                                          make_shared<vector<vector<Pt> > >(torque_di_) :
                                          make_shared<vector<vector<Pt> > >(torque_tri_));
  
  nmer_deets.open( outfile );
  
  for( i = 0; i < nmer->size(); i++)
  {
    mol[0] = nmer->operator[](i)[0];
    mol[1] = nmer->operator[](i)[1];
    mol[2] = nmer->operator[](i)[2];
    if (num == 2 ) dist = _sys_->get_pbc_dist_vec(mol[0], mol[1]).norm();
    else
    {
      vector<double> dists(3);
      Pt com = (_sys_->get_centeri(mol[0]) + _sys_->get_centeri(mol[1]) +
                _sys_->get_centeri(mol[1]))*(1/3.);
      dists[0] = _sys_->get_centeri(mol[0]).dist(com);
      //_sysTemp->get_pbc_dist_vec(0, 1).norm();
      dists[1] = _sys_->get_centeri(mol[1]).dist(com);
      //_sysTemp->get_pbc_dist_vec(0, 2).norm();
      dists[2] = _sys_->get_centeri(mol[2]).dist(com);
      //_sysTemp->get_pbc_dist_vec(1, 2).norm();
      //sort( dists.begin(), dists.end());
      dist = dists[0] + dists[1] + dists[2];
    }
    for ( j = 0; j < num; j++)
    {
      en_nrm = 1.0; //energy_approx_[mol[j]];
      fo_nrm = Pt( 1, 1, 1); //force_approx_[mol[j]];
      print_all[0] = ((abs(en->operator[](i)[j]/en_nrm)<1e-15) ?
                      0 : en->operator[](i)[j]/en_nrm);
      print_all[1] = ((abs(frc->operator[](i)[j].x()/fo_nrm.x())<1e-15) ?
                      0 : frc->operator[](i)[j].x()/fo_nrm.x());
      print_all[2] = ((abs(frc->operator[](i)[j].y()/fo_nrm.y())<1e-15) ?
                      0 : frc->operator[](i)[j].y()/fo_nrm.y());
      print_all[3] = ((abs(frc->operator[](i)[j].z()/fo_nrm.z())<1e-15) ?
                      0 : frc->operator[](i)[j].z()/fo_nrm.z());
      nmer_deets << dist << "\t" << print_all[0] << "\t" << print_all[1];
      nmer_deets << ", " << print_all[2] << ", "<< print_all[3] << endl;
    }
  }
  
  nmer_deets.close();
}


// Two body approximation
void ThreeBody::calcTwoBDEnForTor( )
{
  int i, j, mol;
  
  for( mol = 0; mol < N_; mol++)
  {
    for( i = 0; i < dimer_.size(); i++)
    {
      for( j = 0; j < 2; j++)
      {
        if ( dimer_[i][j] == mol )
        {
          energy_approx_[mol] += energy_di_[i][j];
          force_approx_[mol]   = force_approx_[mol] + force_di_[i][j];
          torque_approx_[mol]  = torque_approx_[mol] + torque_di_[i][j];
        }
      }
    }
  }
  
}

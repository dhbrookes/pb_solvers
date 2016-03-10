//
//  ThreeBody.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/9/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "ThreeBody.h"

ThreeBody::ThreeBody()
{
  
  dimer_.reserve(N_*N_);
  trimer_.reserve(N_*N_*N_);
  
  generatePairsTrips();
}

void ThreeBody::generatePairsTrips()
{
  int i, j, k;
  Pt dist1, dist2, dist3;
  vector<int> temp2(2), temp3(3);
  int M2 = N_*(N_-1)/2; int M3 = N_*(N_-1)*(N_-2)/6;
  
  for( i = 0; i < N_; i++)
  {
    for( j = i+1; j < N_; j++)
    {
      dist1 = _sys_->get_pbc_dist_vec( i, j);
      if (cutoffTBD_ > dist1.norm())
      {
        temp2[0] = i; temp2[1] = j;
        dimer_.push_back(temp2);
      }
      for( k = j + 1; k < N_; k++)
      {
        dist2 = _sys_->get_pbc_dist_vec( i, k);
        dist3 = _sys_->get_pbc_dist_vec( j, k);
        if (cutoffTBD_*2.0 > (dist1.norm() + dist2.norm() + dist3.norm() ))
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
  for ( i = 0; i < dimer_.size(); i++)
  {
    energy_di_[i].resize(2);
    force_di_[i].resize(2);
  }
  
  energy_tri_.resize( trimer_.size());
  force_tri_.resize( trimer_.size());
  for ( i = 0; i < trimer_.size(); i++)
  {
    energy_tri_[i].resize(3);
    force_tri_[i].resize(3);
  }
  
  cout << "Cutoffs implemented, using " << cutoffTBD_ << endl;
  cout << "Max # of di: "  << M2 << " Act used: " << dimer_.size() <<endl;
  cout << "Max # of tri: " << M3 << " Act used: " << trimer_.size() <<endl;
} //end cutoffTBD


// Three body approximation computation
void ThreeBody::solveNmer( int num )
{
  int i, j;
  vector< vector<int> > nmer = ( num == 2) ? dimer_ : trimer_;
  vector< Molecule > mol_temp;
  
  EnergyCalc EnTest;
  ForceCalc FoTest;
  TorqueCalc TorTest;
  System sys_temp;
  ASolver asolv_temp;
  
  for( i = 0; i < nmer.size(); i++)
  {
    cout << "This nmer is: ";
    mol_temp.clear();
    sys_temp = _sys_->get_subsystem(nmer[i]);
    
    asolv_temp = ASolver( num, p_, *_besselCalc_, *_shCalc_, sys_temp);
    asolv_temp.solve_A(1E-5); asolv_temp.solve_gradA(1E-5);
    
    EnTest = EnergyCalc( asolv_temp.get_A(), asolv_temp.calc_L(),
                      _sys_->get_consts(), num, p_);
    FoTest = ForceCalc ( asolv_temp.get_A(), asolv_temp.get_gradA(),
                     asolv_temp.calc_L(), asolv_temp.calc_gradL(),
                     _sys_->get_consts(), num, p_);
    TorTest = TorqueCalc ( *_shCalc_, *_besselCalc_, asolv_temp.calc_gradL(),
                       asolv_temp.get_gamma(), sys_temp.get_consts(),
                       sys_temp, p_);
    
    cout << "Pot, force for all: ";
    for ( j = 0; j < num; j++)
    {
      cout << EnTest.get_omega_i_int(j) << "\t" << FoTest.get_fi(i)[0]
           << ", " << FoTest.get_fi(i)[1]<< ", " << FoTest.get_fi(i)[2]<< "\t";
      if ( num == 2 )
      {
        energy_di_[i][j] = EnTest.get_omega_i_int(j);
        force_di_[i][j] = FoTest.get_fi(i);
      } else {
        energy_tri_[i][j] = EnTest.get_omega_i_int(j);
        force_tri_[i][j] = FoTest.get_fi(i);
      }
    }
    cout << endl;
  }
  cout << num << "mers done " << endl;
  
}
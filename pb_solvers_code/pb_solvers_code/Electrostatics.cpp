//
//  Electrostatics.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/3/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "Electrostatics.h"

Electrostatic::Electrostatic(VecOfMats<cmplx>::type A, System sys,
                             SHCalc shCalc, BesselCalc bCalc, int p, int npts)
: p_(p)
{
  range_min_.resize(3);
  range_max_.resize(3);
  npts_.resize(3);
  step_.resize(3);
  
  A_ = A;
  _sys_ = make_shared<System>(sys);
  _shCalc_ = make_shared<SHCalc> (shCalc);
  _bCalc_  = make_shared<BesselCalc> (bCalc);
  
  for (int i = 0; i < 3; i++)
    npts_[i] = npts;
  
  find_range();
  find_bins();
  
  cout << " This is my range " << range_min_[0] <<  ", " <<range_min_[1] <<  ", "<<range_min_[2] <<  " and max " << range_max_[0] <<  ", " <<range_max_[1] <<  ", "<<range_max_[2] << "  bins "  << step_[0] <<  ", " <<step_[1] <<  ", "<<step_[2] << "  bins "  << npts_[0] <<  ", " <<npts_[1] <<  ", "<<npts_[2] << endl;
  
  compute_pot();
}


void Electrostatic::find_range()
{
  int mol, atom, dim;
  Pt center, curAt;
  double rad;
  double point[3] = {0,0,0}; double max[3] = {0,0,0}; double min[3] = {0,0,0};
  
  int Nmol = _sys_->get_n();
  for ( mol = 0; mol < Nmol; mol++)
  {
    center = _sys_->get_centeri(mol);
    rad    = _sys_->get_radi(mol);
    point[0] = center.x(); point[1] = center.y(); point[2] = center.z();
    for (dim = 0; dim < 3; dim++)
    {
      if(point[dim]-rad < min[dim]) min[dim] = point[dim]-rad;
      if(point[dim]+rad > max[dim]) max[dim] = point[dim]+rad;
    }
    
    for ( atom = 0; atom < _sys_->get_Mi(mol); atom++)
    {
      curAt = center + _sys_->get_posij(mol, atom);
      point[0] = curAt.x(); point[1] = curAt.y(); point[2] = curAt.z();
      for (dim = 0; dim < 3; dim++)
      {
        if(point[dim] < min[dim]) min[dim] = point[dim];
        if(point[dim] > max[dim]) max[dim] = point[dim];
      }
    }
  }
  
  for (dim = 0; dim<3; dim++)
  {
    range_min_[dim] = min[dim] - 10.0;
    range_max_[dim] = max[dim] + 10.0;
  }
}


void Electrostatic::find_bins()
{
  int dim, x, y;
  
  for (dim = 0; dim<3; dim++)
    step_[dim] = (range_max_[dim] - range_min_[dim]) / (double) npts_[dim];
  
  esp_.resize(npts_[0]);
  for ( x = 0; x < npts_[1]; x++)
  {
    esp_[x].resize(npts_[1]);
    for ( y = 0; y < npts_[2]; y++)
      esp_[x][y].resize(npts_[2]);
  }
}

void Electrostatic::print_dx( string dxname )
{
  ofstream dx;
  int xct, yct, zct;
  int ct = 0;
  
  dx.open(dxname);
  dx << "# Data from PBAM Rad run \n# My runname is " << dxname << endl;
  dx << "object 1 class gridpositions counts " << npts_[0]
     << " " << npts_[1] << " " << npts_[2] << endl;
  dx << "origin " << range_min_[0] << " " << range_min_[1]
     << " " << range_min_[2] << endl;
  dx << "delta " << step_[0] << " 0.0e+00 0.0e+00" << endl;
  dx << "delta 0.0e00 " << step_[1] << " 0.0e+00" << endl;
  dx << "delta 0.0e00 0.0e+00 " << step_[2] << endl;
  dx << "object 2 class gridconnections counts " << npts_[0]
    << " " << npts_[1] << " " << npts_[2] << endl;
  dx << "object 3 class array type double rank 0 items "
    << npts_[0]*npts_[1]*npts_[2] << " data follows" << endl;
  
  for ( xct=0; xct<npts_[0]; xct++)
  {
    for ( yct=0; yct<npts_[1]; yct++)
    {
      for ( zct=0; zct<npts_[2]; zct++)
      {
        dx << esp_[xct][yct][zct] << "  ";
        ct++;
        if ((ct % 3) == 0) dx << "\n";
      }
    }
  }
  
  dx << endl; dx << "attribute \"dep\" string \"positions\"" << endl;
  dx << "object \"regular positions regular connections\" class field" << endl;
  dx << "component \"positions\" value 1 " << endl;
  dx << "component \"connections\" value 2" << endl;
  dx << "component \"data\" value 3" << endl;
  dx.close();
  
}

void Electrostatic::print_grid(string dim, double value)
{
  
  
}

void Electrostatic::compute_pot()
{
  int mol, xct, yct, zct;
  Pt center, pos;
  bool cont;
  double rad;
  
  int Nmol = _sys_->get_n();
  double e_s = _sys_->get_consts().get_dielectric_water();
  
  for ( xct=0; xct<npts_[0]; xct++)
  {
    cout  << range_min_[0]+xct*step_[0] << " ..  " << endl ;
    for ( yct=0; yct<npts_[1]; yct++)
    {
      for ( zct=0; zct<npts_[2]; zct++)
      {
        cont = false;
        for ( mol = 0; mol < Nmol; mol++)
        {
          center = _sys_->get_centeri(mol);
          rad    = _sys_->get_radi(mol);
          pos = Pt(range_min_[0]+xct*step_[0], range_min_[1]+yct*step_[1],
                   range_min_[2]+zct*step_[2]);
          
          if ((pos - center).norm2() < rad*rad)
          {
            cont = true;
            break;
          }
        }
        
        if (cont)
          esp_[xct][yct][zct] = 0.0;
        else
          esp_[xct][yct][zct] = compute_pot_at( pos.x(), pos.y(), pos.z())/e_s;
      }
    }
  }
}

double Electrostatic::compute_pot_at( double x, double y, double z)
{

  int mol, Nmol      = _sys_->get_n();
  double rad, pot    = 0.0;
  Pt center, dist;
  MyMatrix<cmplx> localK;
  
  for ( mol = 0; mol < Nmol; mol++)
  {
    center = _sys_->get_centeri(mol);
    rad    = _sys_->get_ai(mol);
    dist   = Pt(x,y,z) - center;
    localK = get_local_exp(dist);
    pot += lotan_inner_prod( A_[mol], localK, p_);
  }
 
  return pot;
}

MyMatrix<cmplx> Electrostatic::get_local_exp( Pt dist )
{
  int n, m;
  double lambda = _sys_->get_lambda();
  double kap    = _sys_->get_consts().get_kappa();
  double expKR;
  vector<double> bessK;
  MyMatrix<cmplx> localK(p_, 2*p_);
  
  bessK = _bCalc_->calc_mbfK(p_, kap*dist.r());
  expKR = exp( - kap * dist.r()) / dist.r();
  _shCalc_->calc_sh(dist.theta(),dist.phi());
  
  for ( n = 0; n < p_; n++)
  {
    for ( m = -n; m <= n; m++)
    {
      localK.set_val( n, m+p_, (pow( lambda/dist.r(), n) * expKR *
                                _shCalc_->get_result(n, m) * bessK[n]));
    }
  }
  
  return localK;
}

double Electrostatic::lotan_inner_prod(MyMatrix<cmplx> U, MyMatrix<cmplx> V,
                                       int p)
{
  double ip = 0;
  int n, m, mT;
  for (n = 0; n < p; n++)
  {
    for (m = -n; m <= n; m++)
    {
      mT = (m < 0) ? -1*m : m;
      ip += U(n, mT+p_).real()*V(n, mT+p_).real()
      + U(n, mT+p_).imag()*V(n, mT+p_).imag();
    }
  }
  return ip;
}


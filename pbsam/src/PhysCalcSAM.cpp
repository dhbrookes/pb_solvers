//
//  PhysCalc.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "PhysCalcSAM.h"


double EnergyCalcSAM::calc_energy(shared_ptr<HMatrix> H,
                                  shared_ptr<LHNMatrix> LHN)
{
  double E(0), eint(0);
  for (int k = 0; k < H->get_ns(); k++)
  {
//    cout << "This is lhn" << endl;
//    LHN->print_kmat(k);
//
//    cout << "This is h" << endl;
//    H->print_kmat(k);

    eint = 0;
    for (int n = 0; n < H->get_p(); n++)
      for (int m = -n; m < n+1; m++)
        eint += (LHN->get_mat_knm(k, n, m).real()*H->get_mat_knm(k, n, m).real()
              +LHN->get_mat_knm(k, n, m).imag()*H->get_mat_knm(k, n, m).imag());
      E += eint;
  }
  return E*epsS;
}

void EnergyCalcSAM::calc_all_energy(vector<shared_ptr<HMatrix> > H,
                                    vector<shared_ptr<LHNMatrix> > LHN)
{
  for (int i = 0; i < H.size(); i++)
  {
    (*omega_)[i] = calc_energy(H[i], LHN[i]);
  //cout << "Molecule energy " << i << " : " << get_ei(i) << endl;
  }
}

void ForceCalcSAM::calc_fI(int I, shared_ptr<HMatrix> H,
                        shared_ptr<LHNMatrix> LHN,
                        shared_ptr<GradHMatrix> dH,
                        shared_ptr<GradLHNMatrix> dLHN)
{
  Pt fIk, tot, inner1, inner2;
  for (int k = 0; k < H->get_ns(); k++)
  {
    forces_[I][k] = 0.0;
    for (int n = 0; n < H->get_p(); n++)
      for (int m = -n; m < n+1; m++)
      {
        for (int d = 0; d < 3; d++)
        {
          inner1.set_cart(d, dLHN->get_mat_knm(k, n, m).get_cart(d).real()
                          * H->get_mat_knm(k, n, m).real() +
                          dLHN->get_mat_knm(k, n, m).get_cart(d).imag()
                          * H->get_mat_knm(k, n, m).imag() );
          inner2.set_cart(d, dH->get_mat_knm(k, n, m).get_cart(d).real()
                          * LHN->get_mat_knm(k, n, m).real() +
                          dH->get_mat_knm(k, n, m).get_cart(d).imag()
                          * LHN->get_mat_knm(k, n, m).imag());
        }
        fIk = inner1 + inner2;
        fIk *= -eps_s_;
        forces_[I][k] += fIk;
      }
    tot += forces_[I][k];
  }
  (*_F_)[I] = tot;
}

void ForceCalcSAM::calc_all_f(vector<shared_ptr<HMatrix> > H,
                           vector<shared_ptr<LHNMatrix> > LHN,
                           vector<vector<shared_ptr<GradHMatrix> > > dH,
                           vector<vector<shared_ptr<GradLHNMatrix> > > dLHN)
{
  for (int i = 0; i < H.size(); i++)
  {
    forces_[i].resize(ks_[i]);
    calc_fI(i, H[i], LHN[i], dH[i][i], dLHN[i][i]);
  //cout << "This is f i " << i << ": " << get_forcei(i).x() << " " <<
  //get_forcei(i).y() << " " << get_forcei(i).z() << " " << endl;
  }
}

PhysCalcSAM::PhysCalcSAM(shared_ptr<Solver> _solv,
                         shared_ptr<GradSolver> _gradsolv,
                         string outfname, Units unit)
:BasePhysCalc(_solv->get_sys()->get_n(), _solv->get_consts(), outfname, unit),
_solv_(_solv), _gradSolv_(_gradsolv), _sys_(_solv->get_sys())
{
  _sys_ = _solv->get_sys();

  _eCalc_ = make_shared<EnergyCalcSAM>(_sys_->get_n(), _solv->get_consts()->get_dielectric_water());
  _fCalc_ = make_shared<ForceCalcSAM>(_sys_->get_n(), _sys_->get_all_Ik(),
                                   _solv->get_consts()->get_dielectric_water(),
                                   _solv->get_sh(), _solv->get_bessel());
  _torCalc_ = make_shared<TorqueCalcSAM>(_sys_->get_n());
}

void TorqueCalcSAM::calc_all_tau(shared_ptr<SystemSAM> sys,
                              shared_ptr<ForceCalcSAM> fcalc)
{
  for (int i = 0; i < I_; i++)
  {
    (*_tau_)[i] = calc_tauI(i, sys->get_moli(i), fcalc);
  }
}

Pt TorqueCalcSAM::calc_tauI(int i, shared_ptr<BaseMolecule> mol,
                         shared_ptr<ForceCalcSAM> fcalc)
{
  Pt ck, rpk, fp, tauI, tau1;
  vector<int> ch_in_k; // charges in sphere k
  vector<Pt> fk = fcalc->get_all_fIk(i);
  for (int k = 0; k < mol->get_ns(); k++)
  {
    ck = mol->get_centerk(k) - mol->get_cog();
    tau1 = cross_prod(ck, fk[k]);
    tauI += tau1;
  }
  return tauI;
}

Pt TorqueCalcSAM::cross_prod(Pt u, Pt v)
{
  Pt c;
  c.set_x(u[1]*v[2] - u[2]*v[1]);
  c.set_y(u[2]*v[0] - u[0]*v[2]);
  c.set_z(u[0]*v[1] - u[1]*v[0]);
  return c;
}


void PhysCalcSAM::calc_force()
{
  _fCalc_->calc_all_f(_solv_->get_all_H(),
                      _solv_->get_all_LHN(),
                      _gradSolv_->get_gradH_all(),
                      _gradSolv_->get_gradLHN_all());
}

void PhysCalcSAM::calc_energy()
{
  _eCalc_->calc_all_energy(_solv_->get_all_H(), _solv_->get_all_LHN());
}


void PhysCalcSAM::calc_torque()
{
  _torCalc_->calc_all_tau(_sys_, _fCalc_);
}

void PhysCalcSAM::print_all()
{
  int i;
  streambuf * buf;
  ofstream of;

  if(outfname_ != "")
  {
    of.open(outfname_, fstream::in | fstream::out | fstream::app);
    buf = of.rdbuf();
  } else {
    buf = cout.rdbuf();
  }

  ostream out(buf);
  out << "My units are " << unit_ << ". Time: " << _sys_->get_time() << endl;
  for ( i = 0; i < N_; i++)
  {
    Pt mol_pos = _sys_->get_cogi(i);
    out << "Molecule #" << i + 1 << endl;
    out << "\tPOSITION: [" << mol_pos.x() << ", " << mol_pos.y();
    out << ", " << mol_pos.z() << "]" << endl;
    out << "\tENERGY: " << unit_conv_ * get_omegai(i) << endl;

    out << "\tFORCE: " << get_forcei(i).norm() * unit_conv_ << ", [";
    out << get_forcei(i).x() * unit_conv_ << " "
        << get_forcei(i).y() * unit_conv_ << " "
        << get_forcei(i).z() * unit_conv_ << "]"<<endl;
    out << "\tTORQUE: " << get_taui(i).norm() << ", [";
    out << get_taui(i).x() * unit_conv_ << " "
        << get_taui(i).y() * unit_conv_ << " "
        << get_taui(i).z() * unit_conv_ << "]"<<endl;
  }
}

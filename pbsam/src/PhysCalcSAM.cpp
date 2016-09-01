//
//  PhysCalc.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "PhysCalcSAM.h"


double EnergyCalcSAM::calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN)
{
  double E = 0;
  for (int k = 0; k < H->get_ns(); k++)
  {
    for (int n = 0; n < H->get_p(); n++)
      for (int m = -n; m < n+1; m++)
        E += (LHN->get_mat_knm(k, n, m).real()*H->get_mat_knm(k, n, m).real()
              +LHN->get_mat_knm(k, n, m).imag()*H->get_mat_knm(k, n, m).imag());
  }
  return E;
}

void EnergyCalcSAM::calc_all_energy(vector<shared_ptr<HMatrix> > H,
                                           vector<shared_ptr<LHNMatrix> > LHN)
{
  for (int i = 0; i < H.size(); i++)
  {
    (*omega_)[i] = calc_energy(H[i], LHN[i]);
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
  (*totForces_)[I] = tot;
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
  }
}

//Ptx ForceCalc::calc_fp(Pt P, shared_ptr<MoleculeSAM> mol,
//                       shared_ptr<HMatrix> H,
//                       shared_ptr<LHNMatrix> LHN,
//                       shared_ptr<GradHMatrix> dH,
//                       shared_ptr<GradLHNMatrix> dLHN)
//{
//  cmplx Hp;
//  Ptx f, fpk, inner1, inner2, dHp;
//  vector<double> besseli;
//
//  shcalc_->calc_sh(P.theta(), P.phi());
//  
//  for (int k = 0; k < H->get_ns(); k++)
//  {
//    besseli = bcalc_->calc_mbfI(H->get_p(), dH->get_kappa()*mol->get_ak(k));
//    for (int n = 0; n < H->get_p(); n++)
//    {
//      for (int m = - n; m < n+1; m++)
//      {
//        Hp = H->make_hb_Ik(k, P, shcalc_, besseli) * shcalc_->get_result(n, m);
//        dHp = dH->calc_dh_P(P, k, besseli, shcalc_);
//        
//        inner1 = dLHN->get_mat_knm(k, n, m) * Hp;
//        inner2 = dHp * LHN->get_mat_knm(k, n, m);
//        fpk = inner1 * (-1) + inner2;
//        fpk *= -1;
//        f += fpk;
//      }
//    }
//  }
//  return f;
//}


void TorqueCalcSAM::calc_all_tau(shared_ptr<SystemSAM> sys,
                              shared_ptr<ForceCalcSAM> fcalc)
{
  for (int i = 0; i < I_; i++)
  {
    (*torques_)[i] = calc_tauI(i, sys->get_moli(i), fcalc);
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

PhysCalcSAM::PhysCalcSAM(shared_ptr<Solver> _solv, shared_ptr<GradSolver> _gradsolv,
                   string outfname, Units unit)
:_solv_(_solv), _gradSolv_(_gradsolv), outfname_(outfname)
{
  _sys_ = _solv->get_sys();
  
  _eCalc_ = make_shared<EnergyCalcSAM>(_sys_->get_n());
  _fCalc_ = make_shared<ForceCalcSAM>(_sys_->get_n(), _sys_->get_all_Ik(),
                                   _solv->get_consts()->get_dielectric_water(),
                                   _solv->get_sh(), _solv->get_bessel());
  _torCalc_ = make_shared<TorqueCalcSAM>(_sys_->get_n());

  compute_units(_solv->get_consts(), unit);
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

void PhysCalcSAM::compute_units( shared_ptr<Constants> cst, Units unit)
{
  if (unit==INTERNAL)
  {
    unit_ = "Internal";
    unit_conv_ = 1.0;
  } else if (unit == KCALMOL)
  {
    unit_  = "kCal/Mol";
    unit_conv_ = cst->convert_int_to_kcal_mol(1.0);
  } else if (unit == JMOL)
  {
    unit_  = "Joules/Mol";
    unit_conv_ = cst->convert_int_to_jmol(1.0);
  } else if (unit == kT)
  {
    unit_  = "kT";
    unit_conv_ = cst->convert_int_to_kT(1.0);
  }
}


//
//  PhysCalc.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "PhysCalc.h"


cmplx EnergyCalc::calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN)
{
  cmplx E = 0;
  for (int k = 0; k < H->get_ns(); k++)
    for (int n = 0; n < H->get_p(); n++)
      for (int m = - n; m < n+1; m++)
        E += LHN->get_mat_knm(k, n, m) * conj(H->get_mat_knm(k, n, m));
  
  return E;
}


Ptx ForceCalc::calc_fI(shared_ptr<HMatrix> H,
                       shared_ptr<LHNMatrix> LHN,
                       shared_ptr<GradHMatrix> dH,
                       shared_ptr<GradLHNMatrix> dLHN)
{
  Ptx f, fIk, inner1, inner2;
  for (int k = 0; k < H->get_ns(); k++)
    for (int n = 0; n < H->get_p(); n++)
      for (int m = - n; m < n+1; m++)
      {
        inner1 = dLHN->get_mat_knm(k, n, m) * H->get_mat_knm(k, n, m);
        inner2 = dH->get_mat_knm(k, n, m) * LHN->get_mat_knm(k, n, m);
        fIk = inner1 * (-1) + inner2;
        fIk *= -1;
        f += fIk;
      }
  return f;
}

Ptx ForceCalc::calc_fp(Pt P, shared_ptr<Molecule> mol,
                       shared_ptr<HMatrix> H,
                       shared_ptr<LHNMatrix> LHN,
                       shared_ptr<GradHMatrix> dH,
                       shared_ptr<GradLHNMatrix> dLHN)
{

  Ptx f, fpk, inner1, inner2, dHp;
  cmplx Hp;
  shcalc_->calc_sh(P.theta(), P.phi());
  vector<double> besseli;
  for (int k = 0; k < H->get_ns(); k++)
  {
    besseli = bcalc_->calc_mbfI(H->get_p(), dH->get_kappa()*mol->get_ak(k));
    for (int n = 0; n < H->get_p(); n++)
    {
      for (int m = - n; m < n+1; m++)
      {
        Hp = H->make_hb_Ik(k, P, shcalc_, besseli) * shcalc_->get_result(n, m);
        dHp = dH->calc_dh_P(P, k, besseli, shcalc_);
        
        inner1 = dLHN->get_mat_knm(k, n, m) * Hp;
        inner2 = dHp * LHN->get_mat_knm(k, n, m);
        fpk = inner1 * (-1) + inner2;
        fpk *= -1;
        f += fpk;
      }
    }
  }
  return f;
}

Ptx TorqueCalc::calc_tauI(shared_ptr<Molecule> mol, shared_ptr<ForceCalc> fcalc,
                          shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN,
                          shared_ptr<GradHMatrix> dH,
                          shared_ptr<GradLHNMatrix> dLHN)
{
  Pt ck, rpk;
  Ptx fk, fp;
  Ptx tauI, tau1, tau2, tau3;
  vector<int> ch_in_k; // charges in sphere k
  for (int k = 0; k < mol->get_ns(); k++)
  {
    ck = mol->get_centerk(k);
    fk = fcalc->calc_fI(H, LHN, dH, dLHN);
    tau1 = cross_prod(ck, fk);
    
    // calculate for every charge on sphere
    tau2 = Ptx();
    ch_in_k = mol->get_ch_allin_k(k);
    for (int p = 0; p < ch_in_k.size(); p++)
    {
      rpk = mol->get_posj(p);
      tau3 = cross_prod(rpk, fcalc->calc_fp(rpk, mol, H, LHN, dH, dLHN));
      tau2 += tau3;
    }
    tauI += tau1;
    tauI += tau2;
  }
  return tauI;
}

Ptx TorqueCalc::cross_prod(Pt u, Ptx v)
{
  Ptx c;
  c.set_x(u[1]*v[2] - u[2]*v[1]);
  c.set_y(u[2]*v[0] - u[0]*v[2]);
  c.set_z(u[0]*v[1] - u[1]*v[0]);
  return c;
}


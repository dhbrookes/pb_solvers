//
//  PhysCalc.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/22/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "PhysCalc.h"


double EnergyCalc::calc_energy(shared_ptr<HMatrix> H, shared_ptr<LHNMatrix> LHN)
{
  double E = 0;
  for (int k = 0; k < H->get_ns(); k++)
  {
//    cout << "This is H " << endl; H->print_kmat(k);
//    cout << "This is LHN " << endl; LHN->print_kmat(k);
    for (int n = 0; n < H->get_p(); n++)
      for (int m = -n; m < n+1; m++)
        E += (LHN->get_mat_knm(k, n, m).real()*H->get_mat_knm(k, n, m).real()
              +LHN->get_mat_knm(k, n, m).imag()*H->get_mat_knm(k, n, m).imag());
//    cout << "This is E " << E << endl;
    
  }
  return E;
}

vector<double> EnergyCalc::calc_all_energy(vector<shared_ptr<HMatrix> > H,
                                           vector<shared_ptr<LHNMatrix> > LHN)
{
  vector<double> E(H.size(), 0.0);
  for (int i = 0; i < H.size(); i++)
  {
//    cout << "This is i " << i << endl;
    E[i] = calc_energy(H[i], LHN[i]);
  }
  return E;
}

void ForceCalc::calc_fI(int I, shared_ptr<HMatrix> H,
                        shared_ptr<LHNMatrix> LHN,
                        shared_ptr<GradHMatrix> dH,
                        shared_ptr<GradLHNMatrix> dLHN)
{
  Pt fIk, inner1, inner2;
  for (int k = 0; k < H->get_ns(); k++)
  {
    forces_[I][k] = 0.0;
//    cout << "This is H" << endl; H->print_kmat(k);
//    cout << "This is dLHN" << endl; dLHN->print_kmat(k);
//    cout << "This is dH" << endl; dH->print_kmat(k);
//    cout << "This is LHN" << endl; LHN->print_kmat(k);
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
//    cout << "For k " << k << " & H*gLHN " << inn1.x() << ", " << inn1.y()
//    << ", " << inn1.z() << endl;
//    cout << "For k " << k << " & gH*gLHN " << inn2.x() << ", " << inn2.y()
//    << ", " << inn2.z() << endl;
//    cout << "{" <<setprecision(9)<< forces_[I][k].x() << "," << forces_[I][k].y() << "," << forces_[I][k].z() << "},"; //<< endl;
  }
}

void ForceCalc::calc_all_f(vector<shared_ptr<HMatrix> > H,
                           vector<shared_ptr<LHNMatrix> > LHN,
                           vector<vector<shared_ptr<GradHMatrix> > > dH,
                           vector<vector<shared_ptr<GradLHNMatrix> > > dLHN)
{
  for (int i = 0; i < H.size(); i++)
  {
//    cout << "This is i " << i << endl;
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


void TorqueCalc::calc_all_tau(shared_ptr<System> sys,
                              shared_ptr<ForceCalc> fcalc)
{
  for (int i = 0; i < I_; i++)
  {
    //    cout << "This is i " << i << endl;
    torques_[i] = calc_tauI(i, sys->get_MoleculeSAM(i), fcalc);
  }
}

Pt TorqueCalc::calc_tauI(int i, shared_ptr<MoleculeSAM> mol,
                         shared_ptr<ForceCalc> fcalc)
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

Pt TorqueCalc::cross_prod(Pt u, Pt v)
{
  Pt c;
  c.set_x(u[1]*v[2] - u[2]*v[1]);
  c.set_y(u[2]*v[0] - u[0]*v[2]);
  c.set_z(u[0]*v[1] - u[1]*v[0]);
  return c;
}


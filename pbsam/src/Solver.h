//
//  Solver.h
//  pbsam_xcode
//
//  Created by David Brookes on 5/17/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#ifndef Solver_h
#define Solver_h

#include <stdio.h>
#include <iostream>
#include <memory>
#include "Gradsolvmat.h"


/*
 Class the uses the above classes to iteratively solve for the F and H matrices
 */
class Solver
{
protected:
  int p_;
  double kappa_;
  int Ns_tot_; // total number of cg sphere
  
  vector<shared_ptr<EMatrix> >      _E_;
  vector<shared_ptr<LEMatrix> >     _LE_;
  vector<shared_ptr<IEMatrix> >     _IE_;
  
  vector<shared_ptr<LFMatrix> >     _LF_;
  vector<shared_ptr<LHMatrix> >     _LH_;
  vector<shared_ptr<LHNMatrix> >    _LHN_;
  vector<shared_ptr<XFMatrix> >     _XF_;
  vector<shared_ptr<XHMatrix> >     _XH_;
  
  vector<shared_ptr<HMatrix> >      _H_;
  vector<shared_ptr<HMatrix> >      _rotH_; // This is for out of molecule?
  vector<shared_ptr<HMatrix> >      _prevH_;
  vector<shared_ptr<HMatrix> >      _outerH_;
  
  vector<shared_ptr<FMatrix> >      _F_;
  
  shared_ptr<TMatrix>               _T_;
  
  shared_ptr<System>                _sys_;
  shared_ptr<SHCalc>                _shCalc_;
  shared_ptr<BesselCalc>            _bCalc_;
  shared_ptr<Constants>             _consts_;
  shared_ptr<ReExpCoeffsConstants>  _reExConsts_;
  shared_ptr<ExpansionConstants>    _expConsts_;
  
  vector<vector<double> >           dev_sph_Ik_;
  
  double                            mu_; // SCF deviation max
  
  // update prevH and outerH and rotH
  void update_rotH(int I, int k);
  void update_outerH(int I, int k);
  void update_prevH(int I, int k);
  void update_prev_all();
  
  void iter_innerH(int I, int k);
  
public:
  Solver(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
         shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
         int p, bool readImat=false, bool readHF=false,
         vector<vector<string> > imats = {{}},
         vector<vector<vector<string> > > expHF = {{{}}});
  
  // run an iteration and return convergence value
  double iter(int t);
  
  // Update LH, LF, XF, XH sequentially
  void step(int t, int I, int k);
  
  // Use h matrices to compute convergence
  double calc_converge_H(int I, int k, bool inner);
  
  void solve(double tol, int maxiter=10000);
  void solve_inner();
  
  void reset_all();
  
  void update_LHN_all();
  
  vector<shared_ptr<HMatrix> > get_all_H() {return _H_;}
  vector<shared_ptr<FMatrix> > get_all_F() {return _F_;}
  MyMatrix<cmplx> getH_ik(int I, int k) {return _H_[I]->get_mat_k(k);}
  MyMatrix<cmplx> getF_ik(int I, int k) {return _F_[I]->get_mat_k(k);}
  cmplx getH_ik_nm(int I, int k, int n, int m)
                  {return _H_[I]->get_mat_knm(k, n, m);}
  cmplx getF_ik_nm(int I, int k, int n, int m)
                  {return _F_[I]->get_mat_knm(k, n, m);}
  
  vector<shared_ptr<LHNMatrix> >  get_all_LHN() {return _LHN_;}
  cmplx getLHN_ik_nm(int I, int k, int n, int m)
                  {return _LHN_[I]->get_mat_knm(k, n, m);}
  
  vector<shared_ptr<IEMatrix> > get_IE()  { return _IE_; }
  
  vector<vector<int> > get_interpol_list()
  {
    vector<vector<int> > ipol(_sys_->get_n());
    for (int i=0; i<_sys_->get_n(); i++) ipol[i]=_LHN_[i]->get_interPol();
    return ipol;
  }
  
  shared_ptr<TMatrix> get_T()              {return _T_;}
  int get_p()                              {return p_;}
  shared_ptr<System> get_sys()             {return _sys_; }
  shared_ptr<Constants> get_consts()       {return _consts_; }
  shared_ptr<SHCalc> get_sh()              {return _shCalc_;}
  shared_ptr<BesselCalc> get_bessel()      {return _bCalc_;}
};


class GradSolver
{
protected:
  int p_;
  double kappa_;
  int Ns_tot_; // measure of number of spheres in the system
  
  vector<shared_ptr<FMatrix> >      _F_;  // converged solutions for these
  vector<shared_ptr<HMatrix> >      _H_;
  vector<shared_ptr<IEMatrix> >     _IE_;

  shared_ptr<TMatrix>               _T_;
  
  // inner index is molecule number, outer index is index of the molecule
  // that this derivative is with respect to
  vector<vector<shared_ptr<GradFMatrix> > >   dF_;
  vector<vector<shared_ptr<GradHMatrix> > >   dH_;
  vector<vector<shared_ptr<GradHMatrix> > >   prev_dH_;
  vector<vector<shared_ptr<GradHMatrix> > >   outer_dH_;
  vector<vector<shared_ptr<GradWFMatrix> > >  dWF_;
  vector<vector<shared_ptr<GradWHMatrix> > >  dWH_;
  vector<vector<shared_ptr<GradLFMatrix> > >  dLF_;
  vector<vector<shared_ptr<GradLHMatrix> > >  dLH_;
  vector<vector<shared_ptr<GradLHNMatrix> > > dLHN_;
  
  vector<vector<shared_ptr<GradCmplxMolMat> > > gradT_A_; // pre-computed part
  
  shared_ptr<System>                _sys_;
  shared_ptr<SHCalc>                _shCalc_;
  shared_ptr<BesselCalc>            _bCalc_;
  shared_ptr<ExpansionConstants>    _expConsts_;
  shared_ptr<Constants>             _consts_;
  
  vector<vector<int> > interpol_; // whether sphere Ik is w/in 10A of other mol
  vector<vector<double > > dev_sph_Ik_; // record of deviations
  
  void iter_inner_gradH(int I, int wrt, int k, vector<double> &besseli,
                        vector<double> &besselk);
  double calc_converge_gradH( int I, int wrt, int k, bool inner);
  
  void step(int t, int I, int wrt, int k, vector<double> &besseli,
            vector<double> &besselk);
  
  // Updating dHs
  void update_prev_gradH(int I, int wrt, int k);
  void update_outer_gradH(int I, int wrt, int k);
  
public:
  GradSolver(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
             shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
             shared_ptr<TMatrix> _T, vector<shared_ptr<FMatrix> > _F,
             vector<shared_ptr<HMatrix> > _H, vector<shared_ptr<IEMatrix> > _IE,
             vector<vector<int> > interpol,
             shared_ptr<ExpansionConstants> _expConst,int p);
  
  void solve(double tol, int maxiter);
  
  void pre_compute_gradT_A();
  
  shared_ptr<GradHMatrix> get_gradH(int I, int wrt) { return dH_[wrt][I];}
  shared_ptr<GradFMatrix> get_gradF(int I, int wrt) { return dF_[wrt][I];}
  
  vector<vector<shared_ptr<GradHMatrix> > > get_gradH_all() { return dH_;}
  vector<vector<shared_ptr<GradFMatrix> > > get_gradF_all() { return dF_;}
  
  vector<vector<shared_ptr<GradLHNMatrix> > > get_gradLHN_all() {return dLHN_;}
  
  cmplx get_gradH_Ik_nm_d(int I, int wrt, int k, int n, int m, int d)
  { return dH_[wrt][I]->get_mat_knm_d(k, n, m, d); }
  cmplx get_gradF_Ik_nm_d(int I, int wrt, int k, int n, int m, int d)
  { return dF_[wrt][I]->get_mat_knm_d(k, n, m, d); }
  
  Ptx get_gradT_A_Ik_nm(int I, int wrt, int k, int n, int m)
  { return gradT_A_[wrt][I]->get_mat_knm(k, n, m); }
  
  double iter(int t, int wrt);
  
};




#endif /* Solver_h */




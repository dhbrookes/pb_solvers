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
#include "Solvmat.h"
#include "Gradsolvmat.h"


/*
 Class the uses the above classes to iteratively solve for the F and H matrices
 */
class Solver
{
protected:
  int p_;
  double kappa_;
  
  vector<shared_ptr<EMatrix> >      _E_;
  vector<shared_ptr<LEMatrix> >     _LE_;
  vector<shared_ptr<IEMatrix> >     _IE_;
  
  vector<shared_ptr<LFMatrix> >     _LF_;
  vector<shared_ptr<LHMatrix> >     _LH_;
  vector<shared_ptr<LHNMatrix> >    _LHN_;
  vector<shared_ptr<XFMatrix> >     _XF_;
  vector<shared_ptr<XHMatrix> >     _XH_;
  
  vector<shared_ptr<HMatrix> >      _H_;
  vector<shared_ptr<HMatrix> >      _prevH_;
  vector<shared_ptr<HMatrix> >      _outerH_;
  
  vector<shared_ptr<FMatrix> >      _F_;
  vector<shared_ptr<FMatrix> >      _prevF_;
  
  shared_ptr<TMatrix>               _T_;
  
  shared_ptr<System>                _sys_;
  shared_ptr<SHCalc>                _shCalc_;
  shared_ptr<BesselCalc>            _bCalc_;
  shared_ptr<Constants>             _consts_;
  shared_ptr<ReExpCoeffsConstants>  _reExConsts_;
  shared_ptr<ExpansionConstants>    _expConsts_;
  
  // update prevH and prevF
  void update_outerH(int I, int k);
  void update_prevH(int I, int k);
  void update_prevF(int I, int k);
  void update_prev_all();
  
  void update_LHN_all();
  
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
  void step(shared_ptr<Molecule> mol, int I, int k);
  
  // Use h matrices to compute convergence
  double calc_converge_H(int I, int k, bool inner);
  
  void solve(double tol, int maxiter=10000);
  void solve_inner();
  
  void reset_all();
  
  vector<shared_ptr<HMatrix> > get_all_H() {return _H_;}
  vector<shared_ptr<FMatrix> > get_all_F() {return _F_;}
  MyMatrix<cmplx> getH_ik(int I, int k) {return _H_[I]->get_mat_k(k);}
  MyMatrix<cmplx> getF_ik(int I, int k) {return _F_[I]->get_mat_k(k);}
  cmplx getH_ik_nm(int I, int k, int n, int m)
                  {return _H_[I]->get_mat_knm(k, n, m);}
  cmplx getF_ik_nm(int I, int k, int n, int m)
                  {return _F_[I]->get_mat_knm(k, n, m);}
  
};


class GradSolver
{
protected:
  int p_;
  double kappa_;
  
  vector<shared_ptr<FMatrix> >      _F_;  // converged solutions for these
  vector<shared_ptr<HMatrix> >      _H_;
  vector<shared_ptr<IEMatrix> >     _IE_;

  shared_ptr<TMatrix>               _T_;
  
  // inner index is molecule number, outer index is index of the molecule
  // that this derivative is with respect to
  vector<vector<shared_ptr<GradFMatrix> > >   dF_;
  vector<vector<shared_ptr<GradHMatrix> > >   dH_;
  vector<vector<shared_ptr<GradWFMatrix> > >  dWF_;
  vector<vector<shared_ptr<GradWHMatrix> > >  dWH_;
  vector<vector<shared_ptr<GradLFMatrix> > >  dLF_;
  vector<vector<shared_ptr<GradLHMatrix> > >  dLH_;
  vector<vector<shared_ptr<GradLHNMatrix> > > dLHN_;
  
  shared_ptr<System>                _sys_;
  shared_ptr<SHCalc>                _shCalc_;
  shared_ptr<BesselCalc>         _bCalc_;
  shared_ptr<Constants>             _consts_;

public:
  GradSolver(shared_ptr<System> _sys, shared_ptr<Constants> _consts,
             shared_ptr<SHCalc> _shCalc, shared_ptr<BesselCalc> _bCalc,
             shared_ptr<TMatrix> _T, vector<shared_ptr<FMatrix> > _F,
             vector<shared_ptr<HMatrix> > _H, int p);
  
  void solve();
  
};




#endif /* Solver_h */




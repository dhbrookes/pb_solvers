//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"

ASolver::ASolver(const int N, const int p, BesselCalc bcalc,
                 SHCalc shCalc, System sys)
:p_(p), a_avg_(sys.get_lambda()), gamma_(N),delta_(N), E_(N), N_(sys.get_n()),
T_ (N_, N_), A_(N), gradA_(N, N), prevA_(N), prevGradA_(N, N), gradT_A_(N, N),
reExpConsts_(sys.get_consts().get_kappa(), sys.get_lambda(), p),
solvedA_(false)
{
  _sys_ = make_shared<System> (sys);
  _besselCalc_ = make_shared<BesselCalc> (bcalc);
  _shCalc_ = make_shared<SHCalc> (shCalc);
  
  // precompute all SH:
  pre_compute_all_sh();
  
  //precompute gamma, delta, E and T:
  compute_T();
  compute_gamma();
  compute_delta();
  compute_E();
  
  init_A();
  init_gradA();
}

// perform many iterations of the solution for A
void ASolver::solve_A(double prec)
{
  int MAX_POL_ROUNDS = 500;

  double scale_dev = (double)(p_*(p_+1)*0.5);
  double cng = scale_dev;
  int ct = 0;
  
  while((cng/scale_dev) > prec)
  {
    iter();
    cng = calc_change();
    if (ct > MAX_POL_ROUNDS*N_)
    {
      cout << "Polz doesn't converge! dev="<< cng << " " << ct << endl;
      exit(0);
    }
    ct++;
  }
  solvedA_ = true;
}

void ASolver::solve_gradA(double prec)
{
  assert(solvedA_); // must solve a before this
  double scale_dev = (double)(p_*(p_+1)*0.5*3.0);
  double cng;
  int j, ct, MAX_POL_ROUNDS = 500;
  
  pre_compute_gradT_A();

  for ( j = 0; j < N_; j++ )
  {
    ct = 0;
    cng = scale_dev;
    while((cng/scale_dev) > prec)
    {
      grad_iter(j);
      cng = calc_grad_change(j);
      if (ct > MAX_POL_ROUNDS*N_)
      {
        cout << "Polz doesn't converge! dev="<< cng << " " << ct << endl;
        exit(0);
      }
      ct++;
    }
  }
  
  for ( int i = 0; i < N_; i++ )
  {
    for ( j = 0; j < N_; j++ )
    {
      print_dAi(i, j, 5);
    }
  }
  
}

// one iteration of numerical solution for A (eq 51 in Lotan 2006)
void ASolver::iter()
{
  int i, j;
  MyMatrix<cmplx> Z, zj, ai;
  prevA_ = A_;
  bool prev = true;  // want to re-expand previous
  for (i = 0; i <  N_; i++)
  {
    // relevant re-expansions:
    Z = MyMatrix<cmplx> (p_, 2*p_ + 1);
    for (j = 0; j < N_; j++)
    {
      if (i == j) continue;
      zj = re_expandA(i, j, prev);
      Z += zj;
    }
    
    ai = delta_[i] * Z;
    ai += E_[i];
    ai = gamma_[i] * ai;
    A_.set_val(i, ai);
  }
}

// one iteration of numerical solution for grad(A) WRT j (eq 53 in Lotan 2006)
void ASolver::grad_iter(int j)
{
  // Solving for grad_j(A^(i)) by iterating through T^(i,k)
  int i, k, dim;

  // relevant re-expansions (g prefix means gradient):
  VecOfMats<cmplx>::type aij, add;
  MyMatrix<cmplx> gamma_delta;
  bool prev = true; //want to re-expand previous
  for (i = 0; i < N_; i++) // molecule of interest
  {
    prevGradA_ = gradA_;
    aij = VecOfMats<cmplx>::type (3);
    for (dim = 0; dim < 3; dim++)
      aij.set_val(dim, MyMatrix<cmplx>(p_,2*p_+1));
    
    aij = get_gradT_Aij( j, i);
    
    cout << "DTA"  << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[0]( n, m+p_).real();
        double im = aij[0]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[1]( n, m+p_).real();
        double im = aij[1]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[2]( n, m+p_).real();
        double im = aij[2]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    
    for (k = 0; k < N_; k++) // other molecules
    {
      if (k == i) continue;
      cout << "This is i " <<i<<" and to J "<<j<<" and k " << k << endl;
      print_dAi( k, j, 5);
      
      add = re_expand_gradA(i, k, j, prev); // T^(i,k) * grad_j A^(k)
      aij += add;
      
      cout << "After Xform"  << endl;
      for (int n = 0; n < 5; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          double  r = add[0]( n, m+p_).real();
          double im = add[0]( n, m+p_).imag();
          r  = fabs( r) > 1e-15 ?  r : 0;
          im = fabs(im) > 1e-15 ? im : 0;
          cout << "(" << setprecision (9) << r << "," << im << ")  ";
          //      cout << "," << setprecision (9) << r ;
        }
        cout << endl;
      }
      cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
      for (int n = 0; n < 5; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          double  r = add[1]( n, m+p_).real();
          double im = add[1]( n, m+p_).imag();
          r  = fabs( r) > 1e-15 ?  r : 0;
          im = fabs(im) > 1e-15 ? im : 0;
          cout << "(" << setprecision (9) << r << "," << im << ")  ";
          //      cout << "," << setprecision (9) << r ;
        }
        cout << endl;
      }
      cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
      for (int n = 0; n < 5; n++)
      {
        for (int m = 0; m <= n; m++)
        {
          double  r = add[2]( n, m+p_).real();
          double im = add[2]( n, m+p_).imag();
          r  = fabs( r) > 1e-15 ?  r : 0;
          im = fabs(im) > 1e-15 ? im : 0;
          cout << "(" << setprecision (9) << r << "," << im << ")  ";
          //      cout << "," << setprecision (9) << r ;
        }
        cout << endl;
      }
      cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
      
      
    }
    
    cout << "Total aij"  << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[0]( n, m+p_).real();
        double im = aij[0]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[1]( n, m+p_).real();
        double im = aij[1]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for (int n = 0; n < 5; n++)
    {
      for (int m = 0; m <= n; m++)
      {
        double  r = aij[2]( n, m+p_).real();
        double im = aij[2]( n, m+p_).imag();
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "(" << setprecision (9) << r << "," << im << ")  ";
        //      cout << "," << setprecision (9) << r ;
      }
      cout << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~" << endl;
    
    gamma_delta = gamma_[i] * delta_[i];
//    aij[0] *= gamma_delta;
    aij.set_val(0, gamma_delta * aij[0]);
    aij.set_val(1, gamma_delta * aij[1]);
    aij.set_val(2, gamma_delta * aij[2]);
    gradA_.set_val(i, j, aij);
    
    cout << "This is mol "<< i<< " and " << j << " new dA" << endl;
    print_dAi( i, j, 5);
  }
}

double ASolver::calc_grad_change(int wrt)
{
  double change = 0.0;
  int i;
  WhichReEx whichA;
  for (i = 0; i < 3; i++)
  {
    if (i == 0) whichA = DDR;
    else if (i == 1) whichA = DDTHETA;
    else whichA = DDPHI;
    change += calc_change(whichA, wrt);
  }
  return change;
}

double ASolver::calc_change(WhichReEx whichA, int wrt)
{
  int i, k, m;
  double change = 0;
  cmplx prev, curr, a; // intermediate values
  for (i = 0; i < N_; i++)
  {
    for(k = 0; k < p_; k++)
    {
      for(m = 0; m <= k; m++)
      {
        // get value from correct verison of A:
        if (whichA==DDR)
        {
          prev = get_prev_dAdr_ni(i, wrt, k, m);
          curr = get_dAdr_ni(i, wrt, k, m);
        }
        else if (whichA==DDTHETA)
        {
          prev = get_prev_dAdtheta_ni(i, wrt, k, m);
          curr = get_dAdtheta_ni(i, wrt, k, m);
        }
        else if (whichA==DDPHI)
        {
          prev = get_prev_dAdphi_ni(i, wrt, k, m);
          curr = get_dAdphi_ni(i, wrt, k, m);
        }
        else
        {
          prev = get_prevA_ni(i, k, m);
          curr = get_A_ni(i, k, m);
        }
        
        if (prev == cmplx(0.0,0.0) && curr == cmplx(0.0,0.0))
          continue;

        if (fabs(prev.real()) < 1e-30 && fabs(prev.imag()) < 1e-30 )
          a = curr;
        else if (fabs(curr.real()) < 1e-30 && fabs(curr.imag()) < 1e-30)
          a = prev;
        else
          a = 0.5*((prev - curr)/(prev + curr));
        
        change += a.real()*a.real() + a.imag()*a.imag();
      }
    }
  }
  return change;
}

// precompute gradT times A(i,j) for all pairs of molecules
void ASolver::pre_compute_gradT_A()
{
  // Solving for grad_j(A^(i)) by iterating through T^(i,k)
  int i, j, k, dim;
  cmplx sign;
  
  // relevant re-expansions (g prefix means gradient):
  VecOfMats<cmplx>::type gjT_Ai, gTA;
  bool prev = true; //want to re-expand previous
  for (i = 0; i < N_; i++) // molecule of interest
  {
    for (j = 0; j < N_; j++) // gradient of interest
    {
      gjT_Ai = VecOfMats<cmplx>::type (3);
      for (dim = 0; dim < 3; dim++)
        gjT_Ai.set_val(dim, MyMatrix<cmplx>(p_,2*p_+1));
      
      if (j == i)
      {
        for (k = 0; k < N_; k++)
        {
          if ( k == i ) continue;
          gTA  = re_expandA_gradT( i, k, prev); // grad_i T^(i,k) A^(k)
          
          sign = (( k < i ) ? cmplx(-1.0, 0.0) : cmplx(1.0, 0.0));
          for (dim = 0; dim < 3; dim++)
            gTA[dim] = sign*gTA[dim];
        
          gjT_Ai += gTA;
        }
      }
      else if (j > i)
        gjT_Ai = re_expandA_gradT( j, i, prev); // grad_j T^(j,i) A^(i)
      else
      {
        gjT_Ai = re_expandA_gradT( j, i, prev); // grad_j T^(j,i) A^(i)
        for (dim = 0; dim < 3; dim++)
          gjT_Ai[dim] = cmplx(-1.0, 0.0)*gjT_Ai[dim];
      }

      gradT_A_.set_val(i, j, gjT_Ai);
    }
  }
  
  
}

//re-expand element j of A with element (i, j) of T and return results
MyMatrix<cmplx> ASolver::re_expandA(int i, int j, bool prev)
{
  MyMatrix<cmplx> x1, x2, z;
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  x1 = expand_RX(  i, j, whichR, whichA, prev);
  x2 = expand_SX(  i, j, x1, whichS);
  z  = expand_RHX( i, j, x2, whichRH);
  return z;
}


/*
 re-expand element j of grad(A) with element (i, j) of T
 will re-expand the gradient with respect to the wrt input
 */
VecOfMats<cmplx>::type ASolver::re_expand_gradA(int i, int j,
                                                int wrt, bool prev)
{
  MyMatrix<cmplx> x1, x2, z;
  VecOfMats<cmplx>::type Z (3);
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=DDR;

  // first dA/dR
  x1 = expand_RX(  i, j, whichR, whichA, prev, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z  = expand_RHX( i, j, x2, whichRH);
  Z.set_val(0, z);

  // dA/dtheta:
  whichA = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA, prev, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z  = expand_RHX( i, j, x2, whichRH);
  Z.set_val(1, z);
  
  // dA/dphiL
  whichA = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA, prev, wrt);
  x2 = expand_SX(  i, j, x1, whichS);
  z  = expand_RHX( i, j, x2, whichRH);
  Z.set_val(2, z);
  
  return Z;
}

/*
 Re-expand element j of A with element (i, j) of grad(T) and return
 as a 3-element vector containing the results for the re-expansion
 with each element of grad(T)
 */
VecOfMats<cmplx>::type ASolver::re_expandA_gradT(int i, int j, bool prev)
{
  MyMatrix<cmplx> x1, x2, z, z1, z2;
  VecOfMats<cmplx>::type Z (3); // output vector
  WhichReEx whichR=BASE, whichS=BASE, whichRH=BASE, whichA=BASE;
  
  int lowI = i; int hiJ = j;
  if ( i > j ) { lowI = j; hiJ  = i; }
  
  // first find with respect to dT/dr:
  whichS = DDR;
  x1 = expand_RX( i, j, whichR, whichA, prev);
  x2 = expand_SX( i, j, x1, whichS);
  z = expand_RHX( i, j, x2, whichRH);
  Z.set_val(0, z);
  
  // dT/dtheta:
  whichR=BASE;
  whichS = BASE;
  whichRH = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA, prev);
  x2 = expand_SX(  i, j, x1, whichS);
  z1 = expand_RHX( i, j, x2, whichRH);

  whichRH = BASE;
  whichR = DDTHETA;
  x1 = expand_RX(  i, j, whichR, whichA, prev);
  x2 = expand_SX(  i, j, x1, whichS);
  z2 = expand_RHX( i, j, x2, whichRH);
  Z.set_val(1, z1 + z2);
  
  // dT/dphi:
  whichR = BASE;
  whichRH = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA, prev);
  x2 = expand_SX(  i, j, x1, whichS);
  z1 = expand_RHX( i, j, x2, whichRH);
  
  whichRH = BASE;
  whichR = DDPHI;
  x1 = expand_RX(  i, j, whichR, whichA, prev);
  x2 = expand_SX(  i, j, x1, whichS);
  z2 = expand_RHX( i, j, x2, whichRH);
  Z.set_val(2, z1 + z2);

  Z = conv_to_cart(Z, i, j);
  return Z;
}

// perform first part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RX(int i, int j, WhichReEx whichR,
                                   WhichReEx whichA, bool prev, int wrt)
{
  int n, m, s, lowI, hiJ;
  MyMatrix<cmplx> x1(p_, 2*p_ + 1);
  cmplx inter, rval, aval;
  
  lowI = i; hiJ = j;
  if ( i > j ) { lowI = j; hiJ  = i; }
  
//  if ( whichA != BASE )  cout << "In RX"  << endl;
  
  // fill X1:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_(lowI,hiJ).isSingular())
      {
        Pt vec = T_(lowI, hiJ).get_TVec();
        if (whichR == BASE)
        {
          aval = which_aval(whichA, prev, j, n, -m, wrt);
          if (vec.theta() > M_PI/2.0)
            x1.set_val(n, m+p_, (n%2 == 0 ? aval : -aval));
          else x1.set_val(n, m+p_, which_aval(whichA, prev, j, n, m, wrt));
        } else if (whichR == DDTHETA)
        {
          x1 = expand_dRdtheta_sing(i, j, vec.theta(), false);
          return x1;
        }
        else
        {
          x1 = expand_dRdphi_sing(i, j, vec.theta(), false);
          return x1;
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichR == DDPHI)
            rval = T_(lowI, hiJ).get_dr_dphi_val(n, m, s);
          else if (whichR == DDTHETA)
            rval = T_(lowI, hiJ).get_dr_dtheta_val(n, m, s);
          else
            rval = T_(lowI, hiJ).get_rval(n, m, s);

          aval = which_aval(whichA, prev, j, n, s, wrt);
//          if ( ( m == 0) and ( whichA != BASE ))
//          {
//                double  r = aval.real();
//                double im = aval.imag();
//                r  = fabs( r) > 1e-15 ?  r : 0;
//                im = fabs(im) > 1e-15 ? im : 0;
//                cout << "(" << setprecision (9) << r << "," << im << ")  ";
//                //      cout << "," << setprecision (9) << r ;
//            
//          }
          inter += rval * aval;
        } // end s
//        if ( ( m == 0) and ( whichA != BASE )) cout << endl;
        x1.set_val(n, m+p_, inter);
      }
    } // end m
  } //end n
  return x1;
}

// perform second part of T*A and return results
MyMatrix<cmplx> ASolver::expand_SX(int i, int j, MyMatrix<cmplx> x1,
                                   WhichReEx whichS)
{
  int n, m, l, lowI, hiJ;
  cmplx inter, sval;
  MyMatrix<cmplx> x2(p_, 2*p_ + 1);
  
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }

  // fill x2:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      for (l = abs(m); l < p_; l++)
      {
        if ( i < j )
        {
          if (whichS == DDR) sval = T_(lowI, hiJ).get_dsdr_val(n, l, m);
          else sval = T_(lowI, hiJ).get_sval(n, l, m);
        }
        else
        {
          if (whichS == DDR) sval = T_(lowI, hiJ).get_dsdr_val(l, n, m);
          else sval = T_(lowI, hiJ).get_sval(l, n, m);
        }
        inter += sval * x1(l, m+p_);
      } // end l
      x2.set_val(n, m+p_, inter);
    } // end m
  } //end n
  return x2;
}

// perform third part of T*A and return results
MyMatrix<cmplx> ASolver::expand_RHX(int i, int j, MyMatrix<cmplx> x2,
                                    WhichReEx whichRH)
{
  int n, m, s, lowI, hiJ;
  cmplx inter, rval;
  MyMatrix<cmplx> z(p_, 2*p_ + 1);
  
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }

  //fill zj:
  for (n = 0; n < p_; n++)
  {
    for (m = -n; m <= n; m++)
    {
      inter  = 0;
      if (T_(lowI,hiJ).isSingular())
      {
        Pt vec = T_(lowI, hiJ).get_TVec();
        if (whichRH == BASE)
        {
          if (vec.theta() > M_PI/2.0)
            z.set_val(n, m+p_, (n%2 == 0 ? x2(n,-m+p_) : -x2(n,-m+p_)));
          else
            z.set_val(n, m+p_, x2(n,m+p_));
        } else if (whichRH == DDTHETA)
        {
          z = expand_dRdtheta_sing(lowI, hiJ, vec.theta(), x2, true);
          return z;
        }
        else
        {
          z = expand_dRdphi_sing(lowI, hiJ, vec.theta(), x2, true);
          return z;
        }
      } else
      {
        for (s = -n; s <= n; s++)
        {
          if (whichRH == DDPHI)
            rval = T_(lowI, hiJ).get_dr_dphi_val(n, s, m);
          else if (whichRH == DDTHETA)
            rval = T_(lowI, hiJ).get_dr_dtheta_val(n, s, m);
          else
            rval = T_(lowI, hiJ).get_rval(n, s, m);

          inter += conj(rval) * x2(n, s+p_);
        } // end s
        z.set_val(n, m+p_, inter);
      }
    } // end m
  } //end n
  return z;
}

MyMatrix<cmplx> ASolver::expand_dRdtheta_sing(int i,int j,double theta,bool ham)
{
  int lowI, hiJ;
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  return expand_dRdtheta_sing(lowI, hiJ, theta, prevA_[j], ham);
}

MyMatrix<cmplx> ASolver::expand_dRdtheta_sing(int i, int j, double theta,
                                              MyMatrix<cmplx> mat, bool ham)
{
  MyMatrix<cmplx> x(p_, 2*p_ + 1);
  x.set_val( 0, p_, cmplx(0.0, 0.0));
  double rec = (ham ? -1.0 : 1.0);
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      x.set_val( n, p_, rec*cmplx(2.0*T_(i,j).get_prefac_dR_val(n,0,1)*
                              mat(n, 1+p_).real(), 0.0)); // m = 0
      for (int m = 1; m < n; m++)
      {
        x.set_val( n, m+p_, rec*T_(i,j).get_prefac_dR_val(n,m,0)*mat(n,m-1+p_)
                   + rec*T_(i,j).get_prefac_dR_val(n,m,1)*mat(n,m+1+p_));
        x.set_val( n,-m+p_, conj( x( n, m+p_)));
      }
      
      x.set_val(n, n+p_,
                rec*T_(i,j).get_prefac_dR_val( n, n, 0)*mat(n, n-1+p_));
      x.set_val(n,-n+p_, conj( x( n, n+p_)));
    }
  }
  else
  {
    double s = -1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      x.set_val( n, p_, rec*cmplx(2.0*s*T_(i,j).get_prefac_dR_val(n,0,1)*
                              mat(n,1+p_).real(), 0.0)); // m = 0
      for (int m = 1; m < n; m++)
      {
        x.set_val( n, m+p_, rec*s*
                  (T_(i,j).get_prefac_dR_val(n,m,0)*mat(n,-m+1+p_)
                   + T_(i,j).get_prefac_dR_val(n,m,1)*mat(n,-m-1+p_)));
        x.set_val( n,-m+p_, conj( x( n, m+p_)));
    }
    
      x.set_val(n, n+p_,rec*s*T_(i,j).get_prefac_dR_val(n, n,0)*mat(n,-n+1+p_));
      x.set_val(n,-n+p_, conj( x( n, n+p_)));
    }
  }
  return x;
}

MyMatrix<cmplx> ASolver::expand_dRdphi_sing(int i,int j, double theta, bool ham)
{
  int lowI, hiJ;
  lowI = i; hiJ = j;
  if ( i > j )  { lowI = j; hiJ  = i; }
  return expand_dRdphi_sing(lowI, hiJ, theta, prevA_[j], ham);
}

MyMatrix<cmplx> ASolver::expand_dRdphi_sing(int i, int j, double theta,
                                            MyMatrix<cmplx> mat, bool ham)
{
  MyMatrix<cmplx> x(p_, 2*p_ + 1);
  x.set_val( 0, p_, cmplx(0.0, 0.0));
  double rec = ((ham && (theta < M_PI/2)) ? -1.0 : 1.0);
  
  if (theta < M_PI/2)
  {
    for (int n = 1; n < p_; n++)
    {
      x.set_val( n, p_, rec*cmplx(2.0*T_(i,j).get_prefac_dR_val(n,0,1)*
                             mat(n, 1+p_).imag(),0.0));
      for (int m = 1; m < n; m++)
      {
        x.set_val(n, m+p_,
            rec*(cmplx( 0.0, T_(i,j).get_prefac_dR_val(n,m,0))*mat(n,m-1+p_)
                - cmplx( 0.0, T_(i,j).get_prefac_dR_val(n,m,1))*mat(n,m+1+p_)));
        x.set_val( n, -m+p_, conj(x(n, m+p_)));
      }
        
      x.set_val( n, n+p_,
            rec*cmplx( 0.0, T_(i,j).get_prefac_dR_val(n,n,0))*mat(n,n-1+p_));
      x.set_val( n, -n+p_, conj(x( n, n+p_)));
    }
  }
  else
  {
    double s = 1.0;
    for (int n = 1; n < p_; n++, s = -s)
    {
      x.set_val(n, p_, rec*cmplx(2.0*s*T_(i,j).get_prefac_dR_val(n,0,1)*
                             mat(n,1+p_).imag(),0.0));
      for (int m = 1; m < n; m++)
      {
        x.set_val(n, m+p_,
            rec*s*(-cmplx(0.0,T_(i,j).get_prefac_dR_val(n,m,0))*mat(n,-m+1+p_)
            + cmplx(0.0, T_(i,j).get_prefac_dR_val(n,m,1))*mat(n,-m-1+p_)));
        x.set_val( n, -m+p_, conj(x(n, m+p_)));
      }
        
      x.set_val( n,  n+p_, rec*cmplx(0.0,
                          -s*T_(i,j).get_prefac_dR_val(n,n,0))*mat(n,-n+1+p_));
      x.set_val( n, -n+p_, conj(x(n, n+p_)));
    }
  }
  return x;
}

/*
 Return correct aval, given A enum and prev flag, mol index, and nm index
 */
cmplx ASolver::which_aval(WhichReEx whichA, bool prev, int i, int n,
                          int m, int wrt)
{
  cmplx aval;
  if (prev)
  {
    if (whichA == DDR)
      aval = get_prev_dAdr_ni(i, wrt, n, m);
    else if (whichA == DDTHETA)
      aval = get_prev_dAdtheta_ni(i, wrt, n, m);
    else if (whichA == DDPHI)
      aval = get_prev_dAdphi_ni(i, wrt, n, m);
    else
      aval = get_prevA_ni(i, n, m);
  }
  else
  {
    if (whichA == DDR)
      aval = get_dAdr_ni(i, wrt, n, m);
    else if (whichA == DDTHETA)
      aval = get_dAdtheta_ni(i, wrt, n, m);
    else if (whichA == DDPHI)
      aval = get_dAdphi_ni(i, wrt, n, m);
    else
      aval = get_A_ni(i, n, m);
  }
  return aval;
}

/*
 Precompute spherical harmonics for E calculation
 */
void ASolver::pre_compute_all_sh()
{
  all_sh.reserve(N_);
  int i;
  
  for (i = 0; i < N_; i++)
    all_sh.push_back(calc_mol_sh(_sys_->get_molecule(i)));
}

/*
 Calculate the SH matrix for every charge in a molecule
 */
vector<MyMatrix<cmplx> > ASolver::calc_mol_sh(Molecule mol)
{
  vector<MyMatrix<cmplx> > vout;
  vout.reserve(mol.get_m());
  int j;
  double theta, phi;
  Pt pt;
  for (j = 0; j < mol.get_m(); j++)
  {
    pt = mol.get_posj(j);
    theta = pt.theta();
    phi = pt.phi();
    _shCalc_->calc_sh(theta, phi);
    vout.push_back(_shCalc_->get_full_result());
  }
  return vout;
}

/*
 Equation 19--Lotan 2006 page 544
 */
cmplx ASolver::calc_indi_gamma(int i, int n)
{
  double g;  // output
  double kap = _sys_->get_consts().get_kappa();
  double ai = _sys_->get_ai(i);
  double eps_p = _sys_->get_consts().get_dielectric_prot();
  double eps_s = _sys_->get_consts().get_dielectric_water();
  // all bessel function k
  vector<double> bk_all = _besselCalc_->calc_mbfK(n+2, kap*ai);
  double bk1 = bk_all[n];  // bessel k at n
  double bk2 = bk_all[n+1];   // bessel k at n+1
  
  double n_dub = (double) n;
  g = (2.0*n_dub + 1.0) * exp(kap*ai);
  g = g / (((2.0*n_dub + 1.0)*bk2) + (n_dub*bk1*((eps_p / eps_s) - 1.0)));
  return cmplx(g, 0);
}

/*
 Equation 20--Lotan 2006 page 544
 */
cmplx ASolver::calc_indi_delta(int i, int n)
{
  double d;  // output
  double kap = _sys_->get_consts().get_kappa();
  double ai = _sys_->get_ai(i);  // radius
  double eps_p = _sys_->get_consts().get_dielectric_prot();
  double eps_s = _sys_->get_consts().get_dielectric_water();
  // all bessel function I:
  vector<double> bi_all = _besselCalc_->calc_mbfI(n+2, kap*ai);
  double bi1 = bi_all[n];  // bessel i at n
  double bi2 = bi_all[n+1];   // bessel i at n+1
  
  double n_dub = (double) n;
  d = (kap*kap*ai*ai) * (bi2 / (2.0*n_dub + 3.0));
  d += (n_dub * bi1 * (1.0 - (eps_p / eps_s)));
  d *= (ai * pow(ai/a_avg_, 2.0*n_dub)) / (2.0*n_dub+1.0);
  return cmplx(d, 0);
}

/*
 Calculates an E^(i)_(n,m) value
 Equation 13--Lotan 2006 page 543
 */
cmplx ASolver::calc_indi_e(int i, int n, int m)
{
  cmplx e = 0.0;
  int j;
  double q, rho, lambda;
  for (j = 0; j < _sys_->get_Mi(i); j++)
  {
    q = _sys_->get_qij(i, j);
    rho = _sys_->get_posij(i, j).r();
    lambda = pow(_sys_->get_lambda(), n);
    // q_ij * (rho_ij)^n * Y_(n,m)(theta_ij, phi_ij):
    cmplx all_sh_acc = all_sh[i][j](n, abs(m));
    if ( m < 0 )
      all_sh_acc = conj( all_sh_acc );

    e += q * ( pow( rho, n ) / lambda ) * all_sh_acc;
  }
  return e;
}


/*
 Constructs the gamma matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_gamma()
 */
void ASolver::compute_gamma()
{
  MyMatrix<cmplx> gi;
  int i, j;
  for (i = 0; i < N_; i++)
  {
    gi = MyMatrix<cmplx> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      gi.set_val(j, j, calc_indi_gamma(i, j));
    }
    gamma_.set_val(i, gi);
  }
}


/*
 Constructs the delta matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_delta()
 */
void ASolver::compute_delta()
{
  MyMatrix<cmplx> di;
  int i, j;
  
  for (i = 0; i < N_; i++)
  {
    di = MyMatrix<cmplx> (p_, p_);
    for(j = 0; j < p_; j++)
    {
      di.set_val(j, j, calc_indi_delta(i, j));
    }
    delta_.set_val(i, di);
  }
}

/*
 Constructs the E vector, which contains a matrix for each molecule
 that defines the multipole expansion of that molecule. The values
 of the inner matrices are calculated in calc_indi_e()
 */
void ASolver::compute_E()
{
  int i, n, m;
  MyMatrix<cmplx> ei;
  for (i = 0; i < N_; i++)
  {
    // m goes from -n to n so you need 2*p columns:
    ei = MyMatrix<cmplx>(p_, 2*p_ + 1);
    for (n = 0; n < p_; n++)
    {
      for (m = -n; m <= n; m++)
      {
        ei.set_val(n, m + p_, calc_indi_e(i, n, m));
      }
    }
    E_.set_val(i, ei);
  }
}

/*
  Calculate the T matrix (i.e the re expansion coefficients along
  every inter-molecular vector (distance between every molecular center
 */
void ASolver::compute_T()
{
  int i, j;
  Pt v, ci, cj;  // inter molecular vector
  
  for (i = 0; i < N_; i++)
  {
    for (j = i+1; j < N_; j++)
    {
      if (i == j) continue;
      ci = _sys_->get_centeri(i);
      cj = _sys_->get_centeri(j);
      v = ci - cj;
      // calculate spherical harmonics for inter molecular vector:
      double kappa = _sys_->get_consts().get_kappa();
      _shCalc_->calc_sh(v.theta(), v.phi());
      vector<double> besselK = _besselCalc_->calc_mbfK(2*p_, kappa * v.r());
      T_.set_val(i, j, ReExpCoeffs(p_, v, _shCalc_->get_full_result(),
                                   besselK, reExpConsts_,
                                   kappa, _sys_->get_lambda(), true));
    }
  }
}

// Initialize A to Gamma * E
void ASolver::init_A()
{
  int i;
  MyMatrix<cmplx> ai;
  for (i = 0; i < N_; i++)
  {
    // m goes from -n to n so you need 2*p columns:
    ai = MyMatrix<cmplx>(p_, 2*p_ + 1);
    ai = gamma_[i] * E_[i];
    A_.set_val(i, ai);
  }
}

// Initialize grad(A) matrix to the zero matrix
void ASolver::init_gradA()
{
  int i, j, k;
  VecOfMats<cmplx>::type ga_ij;
  MyMatrix<cmplx> wrt;
  for (i = 0; i < N_; i++)
  {
    for (j = 0; j < N_; j++)
    {
      ga_ij = VecOfMats<cmplx>::type (3);
      for (k = 0; k < 3; k++)
      {
        wrt = MyMatrix<cmplx>(p_, 2*p_ + 1);
        ga_ij.set_val(k, wrt);
      }
      gradA_.set_val(i, j, ga_ij);
    }
  }
}


/*
 Calculate L, equation (16) of Lotan 2006
 */
VecOfMats<cmplx>::type ASolver::calc_L()
{
  VecOfMats<cmplx>::type L (N_);
  
  int i, j;
  MyMatrix<cmplx> inner, expand;
  for (i = 0; i < N_; i++)
  {
    inner = MyMatrix<cmplx> (p_, 2*p_ + 1);
    for (j = 0; j < N_; j++)
    {
      if (j == i) continue;
      expand = re_expandA(i, j);
      inner += expand;
    }
    L.set_val(i, inner);
  }
  return L;
}

MyVector<VecOfMats<cmplx>::type > ASolver::calc_gradL()
{
  int i, k;
  MyVector<VecOfMats<cmplx>::type > gradl (N_);
  VecOfMats<cmplx>::type inner1, inner2;
  
  for (i = 0; i < N_; i++) // molecule of interest
  {
    inner1 = get_gradT_Aij( i, i);
    
    for (k = 0; k < N_; k++) // other molecules
    {
      if (k == i) continue;
      inner2 = re_expand_gradA(i, k, i, false); // T^(i,k) * grad_j A^(k)
      inner1 += inner2;
    }
    gradl.set_val(i, inner1);
  }
  return gradl;
}

/*
 Convert derivatives from spherical to cartesian coords
 */
VecOfMats<cmplx>::type ASolver::conv_to_cart( VecOfMats<cmplx>::type dZ,
                                               int i, int j )
{
  int m, n; int lowI = i; int hiJ = j;
  double the, r, phi;
  VecOfMats<cmplx>::type Zcart (3); vector<double> con1(3), con2(3), con3(3);

  if ( i > j ) { lowI = j; hiJ  = i; }
  
  Pt v = T_(lowI, hiJ).get_TVec();
  the = v.theta(); r = v.r(); phi = v.phi();

  if (T_(lowI, hiJ).isSingular())
  {
    double cost = (the < M_PI/2) ? 1.0 : -1.0;
    con1 = {  0.0, cost/r,   0.0};
    con2 = {  0.0,    0.0, 1.0/r};
    con3 = { cost,    0.0,   0.0};
  }else
  {
    con1 = {sin(the)*cos(phi),  cos(the)*cos(phi)/r, -sin(phi)/sin(the)/r};
    con2 = {sin(the)*sin(phi),  cos(the)*sin(phi)/r,  cos(phi)/sin(the)/r};
    con3 = {         cos(the),          -sin(the)/r,                  0.0};
  }
  
  for (n = 0; n < 3; n++)
    Zcart[n] = MyMatrix<cmplx> (p_, 2*p_ + 1);
  
  for ( n = 0; n < p_; n++ )
    for ( m = -n; m <= n; m++ )
    {
      (&Zcart[0])->set_val( n, m+p_, con1[0]*dZ[0](n, m+p_)
                        + con1[1]*dZ[1](n, m+p_)+con1[2]*dZ[2](n, m+p_));
      (&Zcart[1])->set_val( n, m+p_, con2[0]*dZ[0](n, m+p_)
                           + con2[1]*dZ[1](n, m+p_)+con2[2]*dZ[2](n, m+p_));
      (&Zcart[2])->set_val( n, m+p_, con3[0]*dZ[0](n, m+p_)
                           + con3[1]*dZ[1](n, m+p_)+con3[2]*dZ[2](n, m+p_));
    }
  return Zcart;
}

/*
 Print function for E
 */
void ASolver::print_Ei( int i, int p)
{
  cout << "This is my E for molecule " << i << " poles " << p <<  endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_E_ni( i, n, m).real();
      double im = get_E_ni( i, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9)  << r << "," << im << ")  ";
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for A
 */
void ASolver::print_Ai( int i, int p)
{
  cout << "This is my A for molecule " << i << " poles " << p <<  endl;
  for (int n = 0; n < 5; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_A_ni( i, n, m).real();
      double im = get_A_ni( i, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << " (" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dx of molecule i wrt j
 */
void ASolver::print_dAidx( int i, int j, int p)
{
  cout << "dAdx for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdx_ni(i, j, n, m).real();
      double im = get_dAdx_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dy of molecule i wrt j
 */
void ASolver::print_dAidy( int i, int j, int p)
{
  cout << "dAdy for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdy_ni(i, j, n, m).real();
      double im = get_dAdy_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA/dz of molecule i wrt j
 */
void ASolver::print_dAidz( int i, int j, int p)
{
  cout << "dAdz for mol " << i << " wrt mol " << j << " pol " << p << endl;
  for (int n = 0; n < p; n++)
  {
    for (int m = 0; m <= n; m++)
    {
      double  r = get_dAdz_ni(i, j, n, m).real();
      double im = get_dAdz_ni(i, j, n, m).imag();
      r  = fabs( r) > 1e-15 ?  r : 0;
      im = fabs(im) > 1e-15 ? im : 0;
      cout << "(" << setprecision (9) << r << "," << im << ")  ";
//      cout << "," << setprecision (9) << r ;
    }
    cout << endl;
  }
  cout << endl;
}

/*
 Print function for dA of molecule i wrt j
 */
void ASolver::print_dAi( int i, int j, int p)
{
  print_dAidx( i, j, p);
  print_dAidy( i, j, p);
  print_dAidz( i, j, p);
}

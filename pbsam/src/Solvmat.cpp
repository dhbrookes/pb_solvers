//
//  Solvmat.cpp
//  pbsam_xcode
//
//  Created by David Brookes on 6/20/16.
//  Copyright © 2016 David Brookes. All rights reserved.
//

#include "Solvmat.h"


vector<Pt> make_uniform_sph_grid(int m_grid, double r)
{
  vector<Pt> grid (m_grid);
  Pt gp;
  grid[0].set_r(r);
  grid[0].set_theta(0.0);
  grid[0].set_phi(0.0);
  double hk;
  for (int k = 1; k < m_grid+1; k++)
  {
    grid[k-1].set_r(r);
    hk = -1.0 + ((2.0 * ( (double)k-1.0)) / ( (double) m_grid-1.0));
    grid[k-1].set_theta(acos(hk));
    
    
    if ( k == 1 || k == m_grid) grid[k-1].set_phi(0);
    else
      grid[k-1].set_phi(grid[k-2].phi() +
                        (3.6/sqrt(m_grid*(1.0-(hk*hk)))));
    
    while(grid[k-1].phi() > 2*M_PI) grid[k-1].set_phi(grid[k-1].phi()-2*M_PI);
    while(grid[k-1].phi() <-2*M_PI) grid[k-1].set_phi(grid[k-1].phi()+2*M_PI);
  }
  return grid;
}

int calc_n_grid_pts(int p, double r)
{
  // not about how these were chosen, taken directly from old stuff
  int npPerOsc(5), nmin(1000), nmax(72000), nFinal;
  int npts = (int) npPerOsc * 0.5 * p;
  if (npts < nmin) npts = nmin;
  if (npts > nmax) npts = nmax;
  
  int minPts = 200;
  double cutoff(1.5), maxrad(10.0);
  const double gradient = npts / ((maxrad-cutoff)*(maxrad-cutoff));
  if (r < cutoff ) nFinal = minPts;
  else
  {
    nFinal = int ( minPts + gradient * (r-cutoff) * (r-cutoff) );
    if (nFinal > npts) nFinal = npts;
  }
  return nFinal;
}

void ExpansionConstants::compute_coeffs()
{
  int ct = 0;
  // constants
  for(int l=0; l < p_; l++)
  {
    set_const1_l( l, (2.0*l+1.0) / (4.0*M_PI));
    set_const2_l( l,  1.0/ (double)(2*l+1));
    for(int m=0; m <= l; m++)
    {
      imatLoc_[ct] = vector<int> {l, m};
      ct++;
    }
  }
}

double inner_prod( MyMatrix<cmplx> & M1, MyMatrix<cmplx> & M2, int p )
{
  double sum = 0;
  
  for (int n = 0; n < p; n++)
    for (int m = -n; m <= n; m++)
    {
      sum += (M1(n, m+p).real() * M2(n, m+p).real()
              + M1(n, m+p).imag() * M2(n, m+p).imag());
    }
  return sum;
}


ExpansionConstants::ExpansionConstants(int p)
: p_(p), expansionConst1_(p), expansionConst2_(p), imatLoc_(p*p)
{
  compute_coeffs();
}


ComplexMoleculeMatrix::ComplexMoleculeMatrix(int I, int ns, int p)
:p_(p), mat_(ns, MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
}

void ComplexMoleculeMatrix::reset_mat(int k)
{
  for (int i = 0; i < mat_[k].get_nrows(); i++)
  {
    for (int j = 0; j < mat_[k].get_ncols(); j++)
    {
      mat_[k].set_val(i, j, cmplx(0, 0));
    }
  }
}

NumericalMatrix::NumericalMatrix(int I, int ns, int p)
:p_(p), mat_(ns), mat_cmplx_(ns, MyMatrix<cmplx> (p, 2*p+1)), I_(I)
{
}

void NumericalMatrix::reset_mat(int k)
{
  //  for (int k = 0; k < mat_.size(); k++)
  //  {
  //    for (int i = 0; i < mat_[k].size(); i++)
  //    {
  //      mat_[k][i] = 0.0;
  //    }
  //  }

  for (int i = 0; i < mat_cmplx_[k].get_nrows(); i++)
  {
    for (int j = 0; j < mat_cmplx_[k].get_ncols(); j++)
    {
      mat_cmplx_[k].set_val(i, j, cmplx(0, 0));
    }
  }
  
}

EMatrix::EMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void EMatrix::calc_vals(shared_ptr<Molecule> mol, shared_ptr<SHCalc> _shcalc,
                        double eps_in)
{
  cmplx val;
  Pt cen;
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol->get_ns(); k++)
  {
    vector<int> allin = mol->get_ch_allin_k(k);
    for (int alpha=0; alpha < allin.size(); alpha++)
    {
      cen = (mol->get_posj_realspace(allin[alpha]) - mol->get_centerk(k));
      _shcalc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol->get_qj(allin[alpha]);
          r_alpha = cen.r();
          a_k = mol->get_ak(k);
          
          val = _shcalc->get_result(n, m);
          val *= q_alpha / eps_in;
          val *= pow(r_alpha / a_k, n);
          val += get_mat_knm(k, n, m);
          set_mat_knm(k, n, m, val);
        }
      }
    }
  }
}

LEMatrix::LEMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void LEMatrix::calc_vals(shared_ptr<Molecule> mol, shared_ptr<SHCalc> _shcalc,
                         double eps_in)
{
  cmplx val;
  Pt cen;
  
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol->get_ns(); k++)
  {
    vector<int> allout = mol->get_ch_allout_k(k);
    for (int alpha=0; alpha < allout.size(); alpha++)
    {
      cen = mol->get_posj_realspace(allout[alpha]) - mol->get_centerk(k);
      _shcalc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol->get_qj(allout[alpha]);
          r_alpha = cen.r();
          a_k = mol->get_ak(k);
          
          val = _shcalc->get_result(n, m);
          val *= q_alpha / eps_in;
          val *= 1 / r_alpha;
          val *= pow(a_k / r_alpha, n);
          val += get_mat_knm(k, n, m);
          set_mat_knm(k, n, m, val);
        }
      }
    }
  }
}


IEMatrix::IEMatrix(int I, shared_ptr<Molecule> _mol,
                   shared_ptr<SHCalc> _shcalc, int p,
                   shared_ptr<ExpansionConstants> _expconst,
                   bool calc_npts, int npts, bool set_mol)
: p_(p), I_(I),
IE_orig_(_mol->get_ns(), vector<double> (p*p*p*p)),
_expConst_(_expconst), calc_pts_(calc_npts), set_mol_(set_mol),
gridPts_(npts), gridPtLocs_(_mol->get_ns()),
grid_exp_(_mol->get_ns()),grid_bur_(_mol->get_ns())
{
  compute_grid_pts(_mol);
}

void IEMatrix::init_from_file(string imatfile, int k )
{
  IMatFile imat(imatfile, p_);
  set_IE_k(k, imat.get_mat());
}


MyMatrix<double> IEMatrix::get_IE_k(int k)
{
  MyMatrix<double> IMat(p_*p_, p_*p_);
  
  int ind(0);
  for (int l = 0; l < p_; l++)  // rows in old matrix
    for (int s = 0; s < 2*l+1; s++)  //columns in old matrix
      for (int n = 0; n < p_; n++)  // rows in new matrix
        for (int m = 0; m < 2*n+1; m++)  // columns in new matrix
        {
          IMat(n*n+m, l*l+s) = get_IE_k_ind(k, ind);
          ind++;
        }
  
  return IMat;
}

void IEMatrix::compute_grid_pts(shared_ptr<Molecule> _mol)
{
  Pt real_grid;
  vector<Pt> grid_loc;
  for (int k = 0; k < _mol->get_ns(); k++)
  {
    vector<int> grid_exp, grid_bur;
    vector<int> neighs = _mol->get_neighj(k);
    if (calc_pts_) gridPts_ = calc_n_grid_pts(p_, _mol->get_ak(k));
    grid_loc = make_uniform_sph_grid(gridPts_, _mol->get_ak(k));
    if (set_mol_)    _mol->set_gridj(k, grid_loc);
    else gridPtLocs_[k] = grid_loc;
    
    // Loop through points to find the exposed ones
    for (int g = 0; g < gridPts_; g++)
    {
      bool buried = false;
      // loop through other spheres in this molecule to find exposed gd pts:
      for (int k2 = 0; k2 < neighs.size(); k2++)
      {
        real_grid = grid_loc[g] + _mol->get_centerk(k);
        // check if sphere is overlapping another
        if (real_grid.dist(_mol->get_centerk(neighs[k2]))
            < _mol->get_ak(neighs[k2]))
        {
          grid_bur.push_back(g);
          buried = true;
          break;
        }
      }
      if (!buried) grid_exp.push_back(g);
    }
    
    if (set_mol_)
    {
      _mol->set_gridexpj(k, grid_exp);
      _mol->set_gridburj(k, grid_bur);
    } else
    {
      grid_bur_[k] = grid_bur;
      grid_exp_[k] = grid_exp;
    }
    
    grid_exp.clear(); grid_bur.clear();
  }
}

vector<MatOfMats<cmplx>::type >
IEMatrix::compute_integral(shared_ptr<Molecule> _mol,
                           shared_ptr<SHCalc> sh_calc,
                           int k)
{
  int min, grid_tot;
  bool bur;
  vector<MatOfMats<cmplx>::type > Ys(2,
                                     MatOfMats<cmplx>::type(p_, p_,
                                                            MyMatrix<cmplx>
                                                            (p_, p_)));
  if (set_mol_)
  {
    gridPtLocs_[k] = _mol->get_gridj(k);
    grid_bur_[k] = _mol->get_gdpt_burj(k);
    grid_exp_[k] = _mol->get_gdpt_expj(k);
  }
  
  grid_tot = (int) (grid_bur_[k].size() + grid_exp_[k].size());
  
  if ( grid_bur_[k].size() < grid_exp_[k].size() )
  {
    bur = true;
    min = (int)grid_bur_[k].size();
  }
  else
  {
    bur = false;
    min = (int)grid_exp_[k].size();
  }
  for(int h=0; h<min; h++)
  {
    double w=0;
    int ind = (bur) ? grid_bur_[k][h] : grid_exp_[k][h];
    Pt gdpt = gridPtLocs_[k][ind];
    // convert the position relative to the center to spherical coordinates.
    sh_calc->calc_sh(gdpt.theta(), gdpt.phi());
    
    // collect sums for (n,m) rows x (l,s) column
    for(int l = 0; l < p_; l++)
      for(int s = 0; s <= l; s++)
      {
        cmplx Yls, Ynm;
        if(s==0) Yls = complex<double> (sh_calc->get_result(l,0));
        else     Yls = complex<double> (sh_calc->get_result(l,s));
        
        for(int n=0; n<=l; n++)
          for(int m=0; m<=n; m++)
          {
            if( n==l && m > s) break;
            if(m==0) Ynm = complex<double> (sh_calc->get_result(n,0));
            else     Ynm = complex<double> (sh_calc->get_result(n,m));
            
            // integrate using the appropriate integration rules
            if(m==0 && s==0 && (n+l)%2==0) //simpson's rule
            {
              if(h > 0 && h < grid_tot)
              {
                if(h % 2 == 0 ) w = 2.0/3.0; // even
                else w = 4.0/3.0; //odd
              } else w = 1.0/3.0;
            } else w = 1.0; // rectangular
            // Yrr and Yri
            Ys[0](l,s).set_val(n, m, Ys[0](l,s)(n,m)+ w * complex<double>
                               (Yls.real()*Ynm.real(), Yls.real()*Ynm.imag()));
            // Yir and Yii
            Ys[1](l,s).set_val(n, m, Ys[1](l,s)(n,m) + w * complex<double>
                               (Yls.imag()*Ynm.real(), Yls.imag()*Ynm.imag()));
          }//nm
      }//ls
    if( h % 10000 == 0 ) cout <<"completed "<<h<<" points "<<endl;
  }//h
  
  double dA = 4*M_PI / (double)(grid_tot-1);
  double mul = (bur) ? -1.0 : 1.0;
  for(int l=0; l<p_; l++)
  {
    double delta = (bur) ? 1.0/_expConst_->get_const1_l(l) : 0.0;
    for(int s=0; s<=l; s++)
      for(int n=0; n<=l; n++)
        for(int m=0; m<=n; m++)
        {
          if( n==l && m > s ) break;
          if( n==l && m==s )
          {
            if(m==0)
            {
              Ys[0](l,s)(n,m).real(delta + mul*dA*Ys[0](l,s)(n,m).real());
              Ys[1](l,s)(n,m).imag(0.0);
            }
            else
            {
              Ys[0](l,s)(n,m).real(0.5*delta + mul*dA*Ys[0](l,s)(n,m).real());
              Ys[1](l,s)(n,m).imag(0.5*delta + mul*dA*Ys[1](l,s)(n,m).imag());
            }
            Ys[0](l,s)(n,m).imag(Ys[0](l,s)(n,m).imag()*mul*dA);
            Ys[1](l,s)(n,m).real(Ys[1](l,s)(n,m).real()*mul*dA);
          }
          else
          {
            Ys[0](l,s)(n,m) *= (mul*dA);
            Ys[1](l,s)(n,m) *= (mul*dA);
          }
        } // end m
  } // end l
  
  for(int l=0; l<p_; l++)
    for(int s=0; s<=l; s++)
      for(int n=0; n<=l; n++)
        for(int m=0; m<=n; m++)
        {
          if( n==l && s < m)
          {
            vector<int> myind = _expConst_->get_imat_loc(m-s-1);
            Ys[0](l,s).set_val(n, m, Ys[0](l,s+1)(myind[0], myind[1]));
            Ys[1](l,s).set_val(n, m, Ys[1](l,s+1)(myind[0], myind[1]));
          }
        }
  
  return Ys;
} // end compute_integral

// Currently using old format because I dont know whats going on!
// TODO: change to understandable
void IEMatrix::populate_mat(vector<MatOfMats<cmplx>::type >  Ys, int k)
{
  int i = 0;
  int sNeg = -1; int mNeg = -1;
  for(int l=0; l<p_; l++)
  {
    double imat;
    for(int n=0; n<p_; n++) // s=0
    {
      double scl = _expConst_->get_const1_l(n);
      for(int m=0; m<=n; m++)
      {
        bool bUpper = ( n<l || (n==l && m==0) );
        imat = scl * (bUpper? Ys[0](l,0)(n,m).real():Ys[0](n,m)(l,0).real());
        IE_orig_[k][i] = imat;
        i++;
        if ( m != 0 )
        {
          imat  = (bUpper ? Ys[0](l,0)(n,m).imag():Ys[1](n,m)(l,0).real());
          IE_orig_[k][i] = (-1.0 * mNeg * scl * imat);
          i++;
        }
      } //m
    }//n
    
    for(int s=1; s<=l; s++) // s>0 and // (l,s)real columns
    {
      for(int n=0; n<p_; n++)
      {
        double scl = _expConst_->get_const1_l(n);
        for(int m=0; m<=n; m++)
        {
          bool bUpper = ( n<l || (n==l && m<=s) );
          imat = 2.*scl*(bUpper? Ys[0](l,s)(n,m).real():Ys[0](n,m)(l,s).real());
          IE_orig_[k][i] = imat;
          i++;
          if ( m != 0 )
          {
            imat  = (bUpper? Ys[0](l,s)(n,m).imag():Ys[1](n,m)(l,s).real());
            IE_orig_[k][i] = (mNeg * -2.0 * scl * imat);
            i++;
          }
        } //m
      }//n
      
      // (l,s)imag columns
      for(int n=0; n<p_; n++)
      {
        double scl = _expConst_->get_const1_l(n);
        for(int m=0; m<=n; m++)
        {
          bool bUpper = ( n<l || (n==l && m<=s) );
          imat  = (bUpper? Ys[1](l,s)(n,m).real():Ys[0](n,m)(l,s).imag());
          IE_orig_[k][i] =  (sNeg * -2.0 * scl * imat);
          i++;
          if ( m != 0 )
          {
            imat  = (bUpper? Ys[1](l,s)(n,m).imag():Ys[1](n,m)(l,s).imag());
            IE_orig_[k][i] = (2.0 * scl * imat);
            i++;
          }
        } //m
      }//n
    }//s
  }//l
}

void IEMatrix::calc_vals(shared_ptr<Molecule> _mol, shared_ptr<SHCalc> _shcalc)
{
  for (int k = 0; k < _mol->get_ns(); k++)
  {
    vector<MatOfMats<cmplx>::type > Ys = compute_integral(_mol, _shcalc, k);
    populate_mat(Ys, k);
  }
}

void IEMatrix::reset_mat()
{
  for (int k = 0; k < IE_orig_.size(); k++)
  {
    for (int ind = 0; ind < p_*p_*p_*p_; ind++)
      IE_orig_[k][ind] = 0.0;
  }
}

LFMatrix::LFMatrix(int I, int ns, int p)
:NumericalMatrix(I, ns, p)
{
}

void LFMatrix::init(shared_ptr<Molecule> mol, shared_ptr<FMatrix> F,
                    shared_ptr<SHCalc> shcalc, shared_ptr<BesselCalc> bcalc,
                    shared_ptr<ExpansionConstants> _expconst)
{
  double rl, im;
  
  for (int k=0; k< mol->get_ns(); k++)
  {
    vector <int> exp_pts = mol->get_gdpt_expj(k);
    mat_[k].resize( (int) exp_pts.size() );
    double dA = 4 * M_PI / (double) mol->get_gridj(k).size();
    
    for (int h = 0; h < exp_pts.size(); h++)
    {
      double val = 0.0;
      Pt q = mol->get_gridjh(k, exp_pts[h]);
      shcalc->calc_sh(q.theta(), q.phi());
      
      for (int n=0; n<p_; n++)
      {
        for (int m = -n; m <= n; m++)
        {
          rl = shcalc->get_result(n,m).real()*F->get_mat_knm(k,n,m).real();
          im = shcalc->get_result(n,m).imag()*F->get_mat_knm(k,n,m).imag();
          val += dA*_expconst->get_const1_l(n) * ( rl + im );
        }
      }
      set_mat_kh(k, h, val);
    }
  }
}

void LFMatrix::calc_vals(shared_ptr<TMatrix> T, shared_ptr<FMatrix> F,
                         shared_ptr<SHCalc> shcalc,
                         shared_ptr<System> sys, int k)
{
  reset_mat(k);
  MyMatrix<cmplx> reex;
  for (int j = 0; j < T->get_nsi(I_); j++)
  {
    if (j==k) continue;
//    print_kmat(j);
    
    if (T->is_analytic(I_, k, I_, j))
    {
      reex = T->re_expandX(F->get_mat_k(j), I_, k, I_, j );
    }
    else
    {
      reex = T->re_expandX_numeric(get_mat(), I_, k, I_, j, 0.0 );
    }
    
//    cout << "For k " << k << " j: " << j << endl;
//    
//    for (int n2 = 0; n2 < p_; n2++)
//    {
//      for (int m2 = 0; m2 < n2+1; m2++)
//      {
//        cout << reex(n2, m2+p_) ;
//      }
//      cout << endl;
//    }
//    
    mat_cmplx_[k] += reex;
  }
}

//TODO: Is this ever needed?
//cmplx LFMatrix::make_fb_Ij(int I, int j, Pt rb,
//                           shared_ptr<FMatrix> F,
//                           shared_ptr<SHCalc> shcalc)
//{
//  cmplx fbj, fb_inner;
//  shcalc->calc_sh(rb.theta(),rb.phi());
//  
//  //create fbj function
//  fbj = 0;
//  for (int n2 = 0; n2 < p_; n2++)
//  {
//    for (int m2 = -n2; m2 < n2+1; m2++)
//    {
//      fb_inner = ((2 * n2 + 1) / (4 * M_PI)) * get_mat_knm(j, n2, m2);
//      fb_inner *= shcalc->get_result(n2, m2);
//      fbj += fb_inner;
//    }
//  }
//  fbj /= pow(rb.r(), 2); // make f_hat into f
//  return fbj;
//}

LHMatrix::LHMatrix(int I, int ns, int p, double kappa)
:NumericalMatrix(I, ns, p), kappa_(kappa)
{
}

void LHMatrix::init(shared_ptr<Molecule> mol, shared_ptr<HMatrix> H,
                    shared_ptr<SHCalc> shcalc, shared_ptr<BesselCalc> bcalc,
                    shared_ptr<ExpansionConstants> _expconst)
{
  double rl, im;
  
  for (int k=0; k< mol->get_ns(); k++)
  {
    vector <int> exp_pts = mol->get_gdpt_expj(k);
    mat_[k].resize( (int) exp_pts.size() );
    double dA = 4 * M_PI / (double) mol->get_gridj(k).size();
    
    for (int h = 0; h < exp_pts.size(); h++)
    {
      double val = 0.0;
      Pt q = mol->get_gridjh(k, exp_pts[h]);
      shcalc->calc_sh(q.theta(), q.phi());
      vector<double> bessI = bcalc->calc_mbfI(p_+1, kappa_*q.r());
      
      for (int n=0; n<p_; n++)
      {
        for (int m = -n; m <= n; m++)
        {
          rl = shcalc->get_result(n,m).real()*H->get_mat_knm(k,n,m).real();
          im = shcalc->get_result(n,m).imag()*H->get_mat_knm(k,n,m).imag();
          val += dA*_expconst->get_const1_l(n)*( rl + im )/bessI[n];
        }
      }
      set_mat_kh(k, h, val);
    }
  }
}

void LHMatrix::calc_vals(shared_ptr<TMatrix> T, shared_ptr<HMatrix> H, int k)
{
  reset_mat(k);
  MyMatrix<cmplx> reex;

  for (int j = 0; j < T->get_nsi(I_); j++)
  {
    if (j==k) continue;

    if (T->is_analytic(I_, k, I_, j))
      reex = T->re_expandX(H->get_mat_k(j), I_, k, I_, j );
    else
      reex = T->re_expandX_numeric(get_mat(), I_, k, I_, j, kappa_ );
    
//        cout << "For k " << k << " j: " << j << endl;
//    
//        for (int n2 = 0; n2 < p_; n2++)
//        {
//          for (int m2 = 0; m2 < n2+1; m2++)
//          {
//            cout << reex(n2, m2+p_) ;
//          }
//          cout << endl;
//        }

    mat_cmplx_[k] += reex;
  }
}

//TODO: Is this ever needed?
//cmplx LHMatrix::make_hb_Ij(int I, int j, Pt rb,
//                           shared_ptr<HMatrix> H,
//                           shared_ptr<SHCalc> shcalc,
//                           vector<double> besseli)
//{
//  //  vector<double> besseli;
//  cmplx hb_inner, hbj;
//  vector<cmplx> h;
//  
//  shcalc->calc_sh(rb.theta(), rb.phi());
//  //  besseli = bcalc->calc_mbfI(p_+1,
//  //                             kappa_*rb.r());
//  hbj = 0;
//  for (int n = 0; n < p_; n++)
//  {
//    for (int m = -n; m < n+1; m++)
//    {
//      hb_inner = ((2*n+1)/ (4*M_PI)) * H->get_mat_knm(j, n, m) ;
//      hb_inner /= besseli[n];
//      hb_inner *= shcalc->get_result(n, m);
//      hbj += hb_inner;
//    }
//  }
//  hbj /= pow(rb.r(), 2);  // make h_hat into h
//  return hbj;
//}

LHNMatrix::LHNMatrix(int I, int ns, int p)
:ComplexMoleculeMatrix(I, ns, p)
{
}

void LHNMatrix::calc_vals(shared_ptr<TMatrix> T,
                          vector<shared_ptr<HMatrix> > H, int k)
{
  reset_mat(k);
  MyMatrix<cmplx> reex;
  for (int J = 0; J < T->get_nmol(); J++)
  {
    if (J == I_) continue;
    for (int l =0 ; l < T->get_nsi(J); l++)
    {
      reex = T->re_expandX(H[J]->get_mat_k(l), I_, k, J, l);
      mat_[k] += reex;
    }
  }

}

XHMatrix::XHMatrix(int I, int ns, int p,
                   shared_ptr<Molecule> mol,
                   shared_ptr<EMatrix> E,
                   shared_ptr<LEMatrix> LE)
: ComplexMoleculeMatrix(I, ns, p), E_LE_mat_(ns, MyMatrix<cmplx> (p, 2*p+1))
{
  double ak;
  
  for (int k = 0; k < ns; k++)
  {
    ak = mol->get_ak(k);
    for (int n = 0; n < p_; n++)
    {
      for (int m = -n; m < n+1; m++)
      {
        E_LE_mat_[k](n, m+p) = (E->get_mat_knm(k, n, m) +
                                ak*LE->get_mat_knm(k, n, m));
      }
    }
  }
}

void XHMatrix::calc_vals(shared_ptr<Molecule> mol, shared_ptr<BesselCalc> bcalc,
                         shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                         shared_ptr<LHNMatrix> LHN, double kappa, int k)
{
  cmplx inner;
  double ak = mol->get_ak(k);
  vector<double> in_k = bcalc->calc_mbfI(p_, kappa * ak);

  for (int n = 0; n < p_; n++)
  {
    for (int m = - n; m <= n; m++)
    {
      inner = E_LE_mat_[k]( n, m+p_);
      inner += ak*LF->get_mat_knm(k, n, m);
      inner -= ak*in_k[n]*(LH->get_mat_knm(k, n, m));
      set_mat_knm(k, n, m, inner);
    }
  }
}


XFMatrix::XFMatrix(int I, int ns, int p, double eps_in, double eps_out,
                   shared_ptr<Molecule> mol, shared_ptr<EMatrix> E,
                   shared_ptr<LEMatrix> LE)
:ComplexMoleculeMatrix(I, ns, p), eps_(eps_in/eps_out),
E_LE_mat_(ns, MyMatrix<cmplx> (p, 2*p+1))
{
  double ak;
  
  for (int k = 0; k < ns; k++)
  {
    ak = mol->get_ak(k);
    for (int n = 0; n < p_; n++)
    {
      for (int m = -n; m <= n; m++)
      {
        E_LE_mat_[k](n, m+p_) = eps_ * ((double)(n+1) * E->get_mat_knm(k, n, m)
                                        - n * ak * LE->get_mat_knm(k, n, m));
      }
    }
  }
}

void XFMatrix::calc_vals(shared_ptr<Molecule> mol, shared_ptr<BesselCalc> bcalc,
                         shared_ptr<LHMatrix> LH, shared_ptr<LFMatrix> LF,
                         shared_ptr<LHNMatrix> LHN, double kappa, int k)
{
  cmplx inner;
  double ak = mol->get_ak(k);
  vector<double> in_k = bcalc->calc_mbfI(p_+1, kappa * ak);
  
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      inner = (pow(kappa * ak, 2.0) * in_k[n+1])/(double)(2*n+3);
      inner += double(n * in_k[n]);
      inner *= ( ak * LH->get_mat_knm(k, n, m));
      inner += E_LE_mat_[k]( n, m+p_);
      inner -= ( n * eps_ * ak * LF->get_mat_knm(k, n, m));
      set_mat_knm(k, n, m, inner);
    }
  }
  
}

HMatrix::HMatrix(int I, int ns, int p, double kappa)
:ComplexMoleculeMatrix(I, ns, p), kappa_(kappa)
{
}


void HMatrix::init_from_exp(string hfilename, int k)
{
  HFFile hfil( hfilename, p_);
  set_mat_k(k, hfil.get_mat());
}


// Initialize H matrix to E with charges mapped to cg (mol.cgCharges_)
void HMatrix::init(shared_ptr<Molecule> mol, shared_ptr<SHCalc> _sh_calc, double eps_in)
{
  cmplx val;
  Pt cen;
  double r_alpha, a_k, q_alpha;
  for (int k=0; k< mol->get_ns(); k++)
  {
    for (int alpha=0; alpha < mol->get_nc_k(k); alpha++)
    {
      cen = mol->get_posj(mol->get_ch_k_alpha(k, alpha));
      _sh_calc->calc_sh(cen.theta(), cen.phi());
      for (int n = 0; n < p_; n++)
      {
        for (int m = -n; m < n+1; m++)
        {
          q_alpha = mol->get_qj(mol->get_ch_k_alpha(k, alpha));
          r_alpha = cen.r();
          a_k = mol->get_ak(k);
          
          val = _sh_calc->get_result(n, m);
          val *= q_alpha / eps_in;
          val *= pow(r_alpha / a_k, n);
          val += get_mat_knm(k, n, m);
          set_mat_knm(k, n, m, val);
        }
      }
    }
  }
}

void HMatrix::calc_vals(shared_ptr<Molecule> mol,
                        shared_ptr<HMatrix> prev,
                        shared_ptr<XHMatrix> XH,
                        shared_ptr<FMatrix> F,
                        shared_ptr<IEMatrix> IE,
                        shared_ptr<BesselCalc> bcalc, int k)
{
  int ct = 0;
  cmplx inner;
  double ak = mol->get_ak(k);
  MyMatrix<double> h_in(p_*p_, 1), h_out(p_*p_, 1);
  vector<double> bessel_i = bcalc->calc_mbfI(p_+1, kappa_*ak);
  vector<double> bessel_k = bcalc->calc_mbfK(p_+1, kappa_*ak);
  
  for (int l = 0; l < p_; l++)  // rows in old matrix
  {
    for (int s = 0; s < l+1; s++)  //columns in old matrix
    {
      inner = (2.0 * l + 1.0) / bessel_i[l];
      inner -= exp(-kappa_*ak) * bessel_k[l];
      inner *= prev->get_mat_knm(k, l, s);
      inner += F->get_mat_knm(k, l, s);
      inner += XH->get_mat_knm(k, l, s);
      
      h_in(ct,0) = inner.real();
      ct++;
      if ( s > 0 )
      {
        h_in(ct,0) = inner.imag();
        ct++;
      }
    }
  }
  
  //TODO: replace with matMul
  h_out = IE->get_IE_k(k) * h_in;  // Matrix vector multiplication
  
  int ctr(0);
  for (int n = 0; n < p_; n++)  // rows in new matrix
  {
    double scl = bessel_i[n] / (double)(2.0*n+1.0);
    for (int m = 0; m < n+1; m++)  // columns in new matrix
    {
      double hRe, hIm;
      hRe = h_out(ctr,0);
      ctr++;
      
      if ( m > 0 )
      {
        hIm = h_out(ctr,0);
        ctr++;
      } else
        hIm = 0.0;
      
      set_mat_knm(k, n, m, scl * complex<double> (hRe, hIm));
      if ( m > 0 ) set_mat_knm(k, n, -m, scl * complex<double> (hRe, -hIm));
    }
  }
}

cmplx HMatrix::make_hb_Ik(int k, Pt rb,
                          shared_ptr<SHCalc> shcalc,
                          vector<double> besseli)
{
  //  vector<double> besseli;
  cmplx hb_inner, hbj;
  vector<cmplx> h;
  
  shcalc->calc_sh(rb.theta(), rb.phi());
  //  besseli = bcalc->calc_mbfI(p_+1,
  //                             kappa_*rb.r());
  hbj = 0;
  for (int n = 0; n < p_; n++)
  {
    for (int m = -n; m < n+1; m++)
    {
      hb_inner = ((2*n+1)/ (4*M_PI)) * get_mat_knm(k, n, m) ;
      hb_inner /= besseli[n];
      hb_inner *= shcalc->get_result(n, m);
      hbj += hb_inner;
    }
  }
  hbj /= pow(rb.r(), 2);  // make h_hat into h
  return hbj;
}



FMatrix::FMatrix(int I, int ns, int p, double kappa)
:ComplexMoleculeMatrix(I, ns, p), kappa_(kappa)
{
}

void FMatrix::init_from_exp(string ffilename, int k)
{
  HFFile ffil( ffilename, p_);
  set_mat_k(k, ffil.get_mat());
}


void FMatrix::calc_vals(shared_ptr<Molecule> mol,
                        shared_ptr<FMatrix> prev,
                        shared_ptr<XFMatrix> XF,
                        shared_ptr<HMatrix> H,
                        shared_ptr<IEMatrix> IE,
                        shared_ptr<BesselCalc> bcalc, int k, double kappa)
{
  int ct = 0;
  cmplx inner;
  double ak = mol->get_ak(k);
  vector<double> bessel_i = bcalc->calc_mbfI(p_+1, kappa*ak);
  vector<double> bessel_k = bcalc->calc_mbfK(p_+1, kappa*ak);
  MyMatrix<double> f_in(p_*p_, 1), f_out(p_*p_, 1);
  
  for (int l = 0; l < p_; l++)
  {
    for (int s = 0; s < l+1; s++)
    {
      inner = double(l * bessel_k[l]) - ((2.0*l+1.0) * bessel_k[l+1]);
      inner *= exp(-kappa*ak) * H->get_mat_knm(k, l, s);
      inner += XF->get_mat_knm(k, l, s);
      inner += (double(2*l + 1 - l * XF->get_eps()) *
                prev->get_mat_knm(k, l, s));
      f_in(ct,0) = inner.real();
      ct++;
      if ( s > 0 )
      {
        f_in(ct,0) = inner.imag();
        ct++;
      }
    }
  }
  
  
  
  //TODO: replace with matMul
  f_out = IE->get_IE_k(k) * f_in;  // Matrix vector multiplication
  
  int ctr(0);
  for (int n = 0; n < p_; n++)  // rows in new matrix
  {
    double scl = 1.0 / (double)(2.0*n+1.0);
    for (int m = 0; m < n+1; m++)  // columns in new matrix
    {
      double fRe, fIm;
      fRe = f_out(ctr,0);
      ctr++;
      
      if ( m > 0 )
      {
        fIm = f_out(ctr,0);
        ctr++;
      } else
        fIm = 0.0;
      
      set_mat_knm(k, n, m, scl * complex<double> (fRe, fIm));
      if ( m > 0 ) set_mat_knm(k, n, -m, scl * complex<double> (fRe, -fIm));
    }
  }
}
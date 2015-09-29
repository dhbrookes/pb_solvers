//
//  ASolver.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "ASolver.h"


ASolver::ASolver(const vector<double>* a, const int N, const int p,
                 const BesselCalc* _bcalc, SHCalc* _shCalc, const System sys)
:p_(p), _besselCalc_(_bcalc), consts_(sys.get_consts()), gamma_(N, N)
,delta_(N, N), E_(N), _shCalc_(_shCalc), sys_(sys), N_(sys.get_n())
{
    // precompute all SH:
    all_sh.reserve(N_);
    int i;
    for (i = 0; i < N_; i++)
    {
        all_sh.push_back(calc_mol_sh(sys.get_molecule(i)));
    }
    
    //precomput gamma, delta and E:
    compute_gamma();
    compute_delta();
    compute_E();
    
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
    ShPt pt;
    for (j = 0; j < mol.get_m(); j++)
    {
        pt = mol.get_sph_posj(j);
        theta = pt.get_theta();
        phi = pt.get_phi();
        _shCalc_->calc_sh(theta, phi);
        vout.push_back(_shCalc_->get_full_result());
    }
    return vout;
}

/*
 Equation 19--Lotan 2006 page 544
 */
const double ASolver::calc_indi_gamma(int i, int n) const
{
    double g;  // output
    double kap = consts_.get_kappa();
    double ai = sys_.get_ai(i);
    double eps_p = consts_.get_dielectric_prot();  // IS THIS WHAT epsilon_p is?? ##########################!!!!!!!!#####################
    double eps_s = consts_.get_dielectric_water();
    vector<double> bk_all = _besselCalc_->calc_mbfK(n+1, kap*ai); // all bessel function k
    double bk2 = bk_all[n];   // bessel k at n+1
    double bk1 = bk_all[n-1];  // bessel k at n
    
    g = (2*n + 1) * exp(kap*ai);
    g = g / (((2*n + 1)*bk2) + (n*bk1*((eps_p / eps_s) - 1)));
    return g;
}

/*
 Equation 20--Lotan 2006 page 544
 */
const double ASolver::calc_indi_delta(int i, int n) const
{
    double d;  // output
    double kap = consts_.get_kappa();
    double ai = sys_.get_ai(i);  // radius
    double eps_p = consts_.get_dielectric_prot();
    double eps_s = consts_.get_dielectric_water();
    vector<double> bi_all = _besselCalc_->calc_mbfK(n+1, kap*ai); // all bessel function I
    double bi2 = bi_all[n];   // bessel i at n+1
    double bi1 = bi_all[n-1];  // bessel i at n
    
//    d = pow(ai, 2*n+1) / (2*n+1);
    d = (kap*kap*ai*ai) * (bi2 / (2*n + 3));
    d += (n * bi1 * (1 - (eps_p / eps_s)));
    d *= pow(ai, 2*n+1) / (2*n+1);
    return d;
}

/*
 Calculates an E^(i)_(n,m) value
 Equation 13--Lotan 2006 page 543
 */
const cmplx ASolver::calc_indi_e(int i, int n, int m)
{
    cmplx e = 0.0;
    int j;
    double q, rho;
    for (j = 0; j < sys_.get_Mi(i); j++)
    {
        q = sys_.get_qij(i, j);
        rho = sys_.get_sph_posij(i, j).get_r();
        e += q * rho * all_sh[i][j](n, m);  // q_ij * rho_ij * Y_(n,m)(theta_ij, phi_ij)
    }
    
    return e;
}


/*
 Constructs the gamma matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_gamma()
 */
void ASolver::compute_gamma()
{
    MyMatrix<double> gi;
    int i, j;
    for (i = 0; i < N_; i++)
    {
        gi = MyMatrix<double> (p_, p_);
        for(j = 0; j < p_; j++)
        {
            gi.set_val(j, j, calc_indi_gamma(i, j));
        }
        gamma_.set_val(i, i, gi);
    }
}


/*
 Constructs the delta matrix, which is a diagonal matrix of diagonal
 matrices that contains values calculated in calc_indi_delta()
 */
void ASolver::compute_delta()
{
    MyMatrix<double> di;
    int i, j;
    for (i = 0; i < N_; i++)
    {
        di = MyMatrix<double> (p_, p_);
        for(j = 0; j < p_; j++)
        {
            di.set_val(j, j, calc_indi_delta(i, j));
        }
        delta_.set_val(i, i, di);
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
        ei = MyMatrix<cmplx>(p_, 2*p_ + 1);  // m goes from -n to n so you need 2*p columns
        for (n = 0; n < p_; n++)
        {
            for (m = -n; m < n; m++)
            {
                ei.set_val(n, m + p_, calc_indi_e(i, n, m));
            }
        }
        E_.set_val(i, ei);
    }
}






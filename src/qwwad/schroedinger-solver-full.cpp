/**
 *  \file     schroedinger-solver-full.cpp
 *  \author   Jonathan Cooper <jdc.tas@gmail.com>
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for full nonparabolic Hamiltonian
 */

#include "schroedinger-solver-full.h"

#include <gsl/gsl_math.h>
#include "constants.h"
#include "linear-algebra.h"

namespace QWWAD
{
using namespace constants;
/**
 * Build matrix 'A' from general eigenproblem
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverFull::SchroedingerSolverFull(const decltype(_m)     &m,
                                               const decltype(_alpha) &alpha,
                                               const decltype(_V)     &V,
                                               const decltype(_z)     &z,
                                               const unsigned int      nst_max) :
    SchroedingerSolver(V,z,nst_max),
    _m(m),
    _alpha(alpha),
    _A(arma::mat(3*z.size(), 3*z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];
    arma::vec B(3*nz*nz);

    // Declare diagonal views
    arma::vec a_elem(nz-1);
    arma::vec b_elem(nz);
    arma::vec c_elem(nz-1);
    arma::vec d_elem(nz-1);
    arma::vec e_elem(nz);
    arma::vec g_elem(nz);

    double const hBar_dz_sq = hBar*hBar/(dz*dz);

    for(unsigned int i=0; i < nz; i++){
        double m_minus;
        double m_plus;
        double alpha_minus;
        double alpha_plus;
        double V_minus;
        double V_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing if neccessary
        if (i==0 or i==nz-1)
        {
            m_minus = m_plus = _m[i];
            alpha_minus = alpha_plus = _alpha[i];
            V_minus = V_plus = _V[i];
        }
        else{
            m_minus = (_m[i]   + _m[i-1])/2;
            m_plus  = (_m[i+1] + _m[i])/2;
            alpha_minus = (_alpha[i]   + _alpha[i-1])/2;
            alpha_plus  = (_alpha[i+1] + _alpha[i])/2;
            V_minus = (_V[i]   + _V[i-1])/2;
            V_plus  = (_V[i+1] + _V[i])/2;
        }

        // Calculate a points
        if(i!=0)
            a_elem(i-1) = -0.5*hBar_dz_sq*\
                                       (1-alpha_plus*V_plus)/(m_minus*alpha_plus*alpha_minus);

        // Calculate b points
        b_elem(i) = 0.5*hBar_dz_sq/(alpha_plus*alpha_minus)*
            ((1-alpha_minus*V_minus)/m_plus + (1-alpha_plus*V_plus)/m_minus)
            + V[i]*(1 - alpha_minus*V_minus - alpha_plus*V_plus +\
                    alpha_plus*alpha_minus*V_plus*V_minus)/(alpha_plus*alpha_minus);

        // Calculate c points
        if(i!=nz-1)
           c_elem(i) = -0.5*hBar_dz_sq*\
                                       (1-alpha_minus*V_minus)/(m_plus*alpha_plus*alpha_minus);

        // Calculate d points
        if(i!=0)
            d_elem(i-1) = -0.5*hBar_dz_sq/(m_minus*alpha_minus);

        // Calculate e points
        e_elem(i) = 0.5*hBar_dz_sq*
            (1/(m_plus*alpha_plus) + 1/(m_minus*alpha_minus))\
            - (1 - alpha_minus*(V_minus+V[i]) - alpha_plus*(V_plus+V[i])\
                    + alpha_plus*alpha_minus*(V_plus*V_minus+V[i]*V_plus+V[i]*V_minus))/(alpha_plus*alpha_minus);

        // Calcualte g points
        g_elem(i) = -1/alpha_plus - 1/alpha_minus + V_plus + V[i]+V_minus;
    }

    // Declare submatrices
    arma::mat A31(nz,nz);
    A31.diag(-1) = a_elem;
    A31.diag()   = b_elem;
    A31.diag(1)  = c_elem;

    // Note that the A32 block is symmetrical so we reuse the d-elements
    arma::mat A32(nz,nz);
    A32.diag(-1) = d_elem;
    A32.diag(0)  = e_elem;
    A32.diag(1)  = d_elem;

    arma::mat A33(nz,nz);
    A33.diag() = g_elem;

    // Insert submatrices into full Hamiltonian matrix
    _A.submat(0,    nz,     nz-1,   2*nz-1).eye(); // A12
    _A.submat(nz,   2*nz,   2*nz-1, 3*nz-1).eye(); // A23
    _A.submat(2*nz, 0,      3*nz-1, nz-1)   = A31;
    _A.submat(2*nz, nz,     3*nz-1, 2*nz-1) = A32;
    _A.submat(2*nz, 2*nz,   3*nz-1, 3*nz-1) = A33;
}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverFull::calculate()
{
    // Find solutions, including all the unwanted "padding" in the eigenvector
    // that comes from the cubic EVP.  See J. Cooper et al., APL 2010
    const auto solutions_tmp = eigen_general(_A, _V.min(), _V.max(), _nst_max);

    // Now chop off the padding from the eigenvector
    const size_t nst = solutions_tmp.size();
    const size_t nz  = _z.size();

    _solutions.clear();
    
    for(unsigned int ist = 0; ist < nst; ist++)
    {
        const auto E   = solutions_tmp[ist].get_E();

        // We just want the first nz elements of the eigenvector
        const auto psi_full = solutions_tmp[ist].psi_array();
        const std::vector<double> psi(psi_full.begin(),
                                      psi_full.begin() + nz);

        _solutions.emplace_back(E, _z, psi);
    }
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

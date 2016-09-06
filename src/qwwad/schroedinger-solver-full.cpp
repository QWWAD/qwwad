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
    _A(arma::vec(9*z.size()*z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];
    arma::vec B(3*nz*nz);

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
            _A[(2*nz)+i+((i-1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)*\
                                       (1-alpha_plus*V_plus)/(m_minus*alpha_plus*alpha_minus);

        // Calculate b points
        _A[(2*nz)+i+(i*3*nz)] = 0.5*gsl_pow_2(hBar/dz)/(alpha_plus*alpha_minus)*
            ((1-alpha_minus*V_minus)/m_plus + (1-alpha_plus*V_plus)/m_minus)
            + V[i]*(1 - alpha_minus*V_minus - alpha_plus*V_plus +\
                    alpha_plus*alpha_minus*V_plus*V_minus)/(alpha_plus*alpha_minus);

        // Calculate c points
        if(i!=nz-1)
            _A[(2*nz)+i+((i+1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)*\
                                       (1-alpha_minus*V_minus)/(m_plus*alpha_plus*alpha_minus);

        // Calculate d points
        if(i!=0)
            _A[(3*nz*nz)+(2*nz)+i+((i-1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)/(m_minus*alpha_minus);

        // Calculate e points
        _A[(3*nz*nz)+(2*nz)+i+(i*3*nz)] = 0.5*gsl_pow_2(hBar/dz)*
            (1/(m_plus*alpha_plus) + 1/(m_minus*alpha_minus))\
            - (1 - alpha_minus*(V_minus+V[i]) - alpha_plus*(V_plus+V[i])\
                    + alpha_plus*alpha_minus*(V_plus*V_minus+V[i]*V_plus+V[i]*V_minus))/(alpha_plus*alpha_minus);

        // Calculate f points
        if(i!=nz-1)
            _A[(3*nz*nz)+(2*nz)+i+((i+1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)/(m_plus*alpha_plus);

        // Calcualte g points
        _A[(6*nz*nz)+(2*nz)+i+(i*3*nz)] = -1/alpha_plus - 1/alpha_minus + V_plus + V[i]+V_minus;

        // Insert identity matrices
        _A[(3*nz*nz)+i+(i*3*nz)] = 1;
        _A[(6*nz*nz)+nz+i+(i*3*nz)] = 1;
    }
}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverFull::calculate()
{
    // Find solutions, including all the unwanted "padding" in the eigenvector
    // that comes from the cubic EVP.  See J. Cooper et al., APL 2010
    const auto solutions_tmp = eigen_general(&_A[0], _V.min(), _V.max(), 3*_z.size(), _nst_max);

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

        _solutions.push_back(Eigenstate(E, _z, psi));
    }
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

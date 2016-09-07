/**
 *  \file     schroedinger-solver-tridiagonal.cpp
 *  \author   Jonathan Cooper <jdc.tas@gmail.com
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions using tridiagonal matrix
 */

#include "schroedinger-solver-tridiagonal.h"
#include <gsl/gsl_math.h>

#include "constants.h"
#include "linear-algebra.h"

namespace QWWAD
{
using namespace constants;
/**
 * Create tridiagonal Hamiltonian
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverTridiag::SchroedingerSolverTridiag(const decltype(_m) &me,
                                                     const decltype(_V) &V,
                                                     const decltype(_z) &z,
                                                     const unsigned int  nst_max) :
    SchroedingerSolver(V,z,nst_max),
    diag(arma::vec(z.size())),
    sub(arma::vec(z.size()-1))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    for(unsigned int i=0; i<nz; i++){
        double m_minus;
        double m_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing
        if(i==0 || i==nz-1){
            m_minus = m_plus = me[i];
        }
        else{
            m_minus = (me[i] + me[i-1])/2;
            m_plus = (me[i+1] + me[i])/2;
        }

        // Calculate a points
        if(i!=nz-1) sub[i] = -gsl_pow_2(hBar/dz)/(2*m_plus);

        // Calculate b points
        diag[i] = 0.5*gsl_pow_2(hBar/dz)*(m_plus+m_minus)/(m_plus*m_minus) + V[i];
    }
}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverTridiag::calculate()
{
    const auto EVP_solutions =
        _E_cutoff_set ?
        eigen_tridiag(diag, sub, _V.min(), _E_cutoff)
        :
        eigen_tridiag(diag, sub, _V.min(), _V.max(), _nst_max);

    _solutions.clear();

    for (auto st : EVP_solutions)
    {
        const auto E   = st.get_E();
        const auto psi = st.psi_array();
        _solutions.push_back(Eigenstate(E, _z, psi));
    }
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

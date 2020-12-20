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
                                                     const arma::vec    &V,
                                                     const arma::vec    &z,
                                                     const unsigned int  nst_max) :
    diag(arma::zeros(z.size())),
    sub(arma::zeros(z.size()-1))
{
    set_V(V);
    set_z(z);
    set_nst_max(nst_max);

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
auto
SchroedingerSolverTridiag::calculate() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;

    const auto z = get_z();

    // Get limits for search
    const double E_min = get_E_search_min();
    const double E_max = get_E_search_max();

    // Set number of states only if energy limits haven't been specified
    // Note that '0' means that we should find all states in range
    const double nst_max = (get_E_min_set() || get_E_max_set()) ? 0 : get_nst_max();

    const auto EVP_solutions = eigen_tridiag(diag, sub, E_min, E_max, nst_max);

    for (auto st : EVP_solutions)
    {
        const auto E   = st.get_E();
        const auto psi = st.psi_array();
        solutions.emplace_back(E, z, psi);
    }

    return solutions;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

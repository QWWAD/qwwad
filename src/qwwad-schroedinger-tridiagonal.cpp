/**
 *  \file     qwwad-schroedinger-tridiagonal.cpp
 *  \author   Jonathan Cooper <jdc.tas@gmail.com
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions using tridiagonal matrix
 */

#include "qwwad-schroedinger-tridiagonal.h"
#include <gsl/gsl_math.h>

#include "qclsim-constants.h"

namespace Leeds
{
using namespace constants;

/**
 * Create tridiagonal Hamiltonian
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverTridiag::SchroedingerSolverTridiag(const std::valarray<double>& me,
                                                     const std::valarray<double>& V,
                                                     const std::valarray<double>& z,
                                                     const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    diag(std::valarray<double>(z.size())),
    sub(std::valarray<double>(z.size()))
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
    if (_E_cutoff_set)
        _solutions = eigen_tridiag(&diag[0], &sub[0], _V.min(), _E_cutoff, _V.size());
    else
        _solutions = eigen_tridiag(&diag[0], &sub[0], _V.min(), _V.max(), _V.size(), _nst_max);
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

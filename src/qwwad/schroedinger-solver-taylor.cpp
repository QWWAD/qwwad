/**
 *  \file   schroedinger-solver-taylor.cpp
 *  \author Jonathan Cooper <jdc.tas@gmail.com
 *  \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief  Implementatation of Schrodinger solver using Taylor expansion
 */

#include "schroedinger-solver-taylor.h"
#include "constants.h"

namespace QWWAD
{
using namespace constants;
/**
 * Create discretised Hamiltonian for system
 *
 * \param[in] me      Effective mass
 * \param[in] alpha   Band nonparabolicity
 * \param[in] V       Confining potential
 * \param[in] z       Spatial coordinates
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverTaylor::SchroedingerSolverTaylor(const std::valarray<double> &me,
                                       const std::valarray<double> &alpha,
                                       const std::valarray<double> &V,
                                       const std::valarray<double> &z,
                                       const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    AB(std::valarray<double>(2*z.size())),
    BB(std::valarray<double>(2*z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    for(unsigned int i=0; i<nz; i++){
        double m_minus;
        double m_plus;
        double alpha_minus;
        double alpha_plus;
        double V_minus;
        double V_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing
        if(i==0 or i==nz-1)
        {
            m_minus     = m_plus     = me[i];
            alpha_minus = alpha_plus = alpha[i];
            V_minus     = V_plus     = V[i];
        }
        else
        {
            m_minus     = (me[i] + me[i-1])/2;
            m_plus      = (me[i+1] + me[i])/2;
            alpha_minus = (alpha[i] + alpha[i-1])/2;
            alpha_plus  = (alpha[i+1] + alpha[i])/2;
            V_minus     = (V[i] + V[i-1])/2;
            V_plus      = (V[i+1] + V[i])/2;
        }

        if(i!=nz-1)
        {
            // Calculate a points
            AB[2*(i+1)] = -0.5*gsl_pow_2(hBar/dz)*(1+alpha_plus*V_plus)/m_plus;

            // Calculate d points
            BB[2*(i+1)] = -0.5*gsl_pow_2(hBar/dz)*alpha_plus/m_plus;
        }

        // Calculate b points
        AB[1+(2*i)] = 0.5*gsl_pow_2(hBar/dz)*((1.0+alpha_plus*V_plus)/m_plus + (1.0+alpha_minus*V_minus)/m_minus) + V[i];

        // Calculate e points
        BB[1+(2*i)] = 0.5*gsl_pow_2(hBar/dz)*(alpha_plus/m_plus + alpha_minus/m_minus) + 1;
    }
}

/**
 * Find solutions to Schroedinger's equation for this Hamiltonian
 */
void SchroedingerSolverTaylor::calculate()
{
    _solutions = eigen_banded(&AB[0], &BB[0], _V.min(), _V.max(), _V.size(), _nst_max);
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

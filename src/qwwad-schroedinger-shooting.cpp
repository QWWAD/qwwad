/**
 *  \file     qwwad-schroedinger-shooting.cpp
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \author   Paul Harrison <p.harrison@shu.ac.uk>
 *  \brief    Implementatation of Schrodinger solver using shooting method
 */
#include "qwwad-schroedinger-shooting.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "qwwad/constants.h"

using namespace QWWAD;
using namespace constants;

namespace Leeds
{
/**
 * \brief Set system parameters for solver
 *
 * \param[in] me      Band-edge effective mass [kg]
 * \param[in] alpha   Nonparabolicity parameter [1/J]
 * \param[in] V       Band-edge potential [J]
 * \param[in] z       Spatial locations [m]
 * \param[in] dE      Minimum energy separation between states [J]
 * \param[in] nst_max Maximum number of states to find
 */
SchroedingerSolverShooting::SchroedingerSolverShooting(const std::valarray<double>& me,
                                                       const std::valarray<double>& alpha,
                                                       const std::valarray<double>& V,
                                                       const std::valarray<double>& z,
                                                       const double                 dE,
                                                       const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    _me(me),
    _alpha(alpha),
    _dE(dE)
{}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverShooting::calculate()
{
    double Elo=_V.min() + _dE;    // first energy estimate

    gsl_function f;
    f.function = &psi_at_inf;
    f.params   = this;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    for(unsigned int ist=0;
            (_nst_max > 0  && ist < _nst_max) || // Continue if max. states is specified & not exceeded
            (_nst_max == 0 && gsl_fcmp(Elo, _V.max(), e/1e12) == -1); // Or if max. states is NOT specified and we're below the max potential
    ++ist)
    {
        // Shift the lower estimate up past the last state we found
        if(ist > 0)
        {
            const double E_last = _solutions[ist-1].get_E();
            Elo = E_last + _dE;
        }

        // Value for y=f(x) at bottom of search range
        const double y1 = GSL_FN_EVAL(&f,Elo);

        // Find the range in which the solution lies by incrementing the
        // upper limit of the search range until the function changes sign.
        // Since this coarse search can require many iterations, we keep the
        // lower limit fixed to minimise the amount of computation at each step.
        // The Brent algorithm is extremely fast, so it really doesn't matter that
        // the range we find here is large.
        //
        // Note the end stop to prevent infinite loop in absence of solution
        //
        // TODO: Make the cut-off configurable
        double y2 = y1;
        double Ehi = Elo;
        do
        {
            Ehi += _dE;
            y2=GSL_FN_EVAL(&f, Ehi);
        }while(y1*y2>0);

        double E = (Elo + Ehi)/2;
        gsl_root_fsolver_set(solver, &f, Elo, Ehi);
        int status = 0;

        // Improve the estimate of the solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            E   = gsl_root_fsolver_root(solver);
            Elo = gsl_root_fsolver_x_lower(solver);
            Ehi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(Elo, Ehi, 1e-12*e, 0);
        }while(status == GSL_CONTINUE);

        // Stop if we've exceeded the cut-off energy
        if(_E_cutoff_set && gsl_fcmp(E, _E_cutoff, e*1e-12) == 1)
            break;

        std::valarray<double> psi(_z.size());
        const double psi_inf = shoot_wavefunction(psi, E);

        _solutions.push_back(State(E,psi));

        // Check that wavefunction is tightly bound
        // TODO: Implement a better check
        if(gsl_fcmp(fabs(psi_inf), 0, 1) == 1)
            throw "Warning: Wavefunction is not tightly bound";
    }
}

/**
 * \brief Find the wavefunction just beyond the right-hand side of the system
 *
 * \details This function returns the value of the wavefunction (psi)
 *          at +infinity for a given value of the energy.  The solution
 *          to the energy occurs for psi(+infinity)=0.
 *
 * \param[in] E      Energy [J]
 * \param[in] params SchroedingerSolverShooting object for which to solve
 *
 * \returns The wavefunction amplitude immediately to the right of the structure
 */
double SchroedingerSolverShooting::psi_at_inf(double  E,
                                              void   *params)
{
    const SchroedingerSolverShooting *se = reinterpret_cast<SchroedingerSolverShooting *>(params);
    std::valarray<double> psi(se->get_z().size());

    const double psi_inf = se->shoot_wavefunction(psi, E);
    return psi_inf;
}

/**
 * \brief Computes wavefunction iteratively from left to right of structure
 *
 * \details The value of the wavefunction is taken to be zero at the point
 *          immediately to the left of the potential profile. Subsequent
 *          values are computed using QWWAD3, Eq. 3.53.
 *
 * \param[out] wf      Array to which wavefunction will be written [m^{-1/2}]
 * \param[in]  E       Energy at which to compute wavefunction
 *
 * \returns The wavefunction amplitude at the point immediately to the right of the structure
 */
double SchroedingerSolverShooting::shoot_wavefunction(std::valarray<double> &wf,
                                                      const double           E) const
{
    const size_t nz = _z.size();
    wf.resize(nz);
    const double dz = _z[1] - _z[0];

    // Recalculate effective mass with non-parabolicity at this energy
    std::valarray<double> m = _me*(1.0+_alpha*(E-_V));

    // boundary conditions (psi[-1] = psi[n] = 0)
    wf[0]   = 1.0;
    double wf_next = 1.0;

    for(unsigned int i=0; i < nz; i++) // last potential not used
    {
        double wf_prev = 0;

        // Compute m(z + dz/2)
        double m_prev = 0.0;
        double m_next = 0.0;

        if(i != 0)
        {
            wf_prev = wf[i-1];
            m_prev = (m[i] + m[i-1])/2.0;
        }
        else
        {
            m_prev = m[i];
        }

        if(i != nz - 1)
            m_next = (m[i] + m[i+1])/2.0;
        else
            m_next = m[i];

        wf_next = (2*m_next*dz*dz/hBar/hBar*(_V[i]-E)+
                1.0 + m_next/m_prev)*wf[i]
                - wf_prev * m_next/m_prev;
        wf_prev += 0;

        if(i != nz-1) wf[i+1] = wf_next;
    } 

    // Normalise the wavefunction
    const std::valarray<double> pd = wf*wf;
    const double norm = trapz(pd,dz);
    wf/= sqrt(norm);

    return wf_next/sqrt(norm);
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

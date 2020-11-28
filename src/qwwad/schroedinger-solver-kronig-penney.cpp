/**
 *  \file   schroedinger-solver-kronig-penney.cpp
 *  \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief  Implementatation of Schrodinger solver functions for Kronig-Penney superlattice
 */

#include "schroedinger-solver-kronig-penney.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "constants.h"
#include "maths-helpers.h"
#include <complex>
#include <armadillo>

namespace QWWAD
{
using namespace constants;
SchroedingerSolverKronigPenney::SchroedingerSolverKronigPenney(const double l_w,
                                                               const double l_b,
                                                               const double V0,
                                                               const double m_w,
                                                               const double m_b,
                                                               const double k,
                                                               const size_t nz,
                                                               const size_t nper,
                                                               const unsigned int nst_max) :
    SchroedingerSolver(arma::zeros(nz*nper),
                       arma::zeros(nz*nper),
                       nst_max),
    _l_w(l_w),
    _l_b(l_b),
    _V0(V0),
    _m_w(m_w),
    _m_b(m_b),
    _k(k)
{
    const auto L  = l_b + l_w;
    const auto dz = L*nper/(nz*nper+1);

    // Generate potential profile (well, then barrier)
    for(unsigned int iper=0; iper < nper; ++iper)
    {
        for (unsigned int iz=0; iz<nz; iz++)
        {
            _z(iz + iper*nz) = iz*dz + iper*L;

            // Fill in barrier potential
            if(_z(iz) > l_w)
                _V(iz + iper*nz) = _V0;
        }
    }
}

/**
 * \brief Finds the right-hand side of the matching equation for a finite well
 *
 * \returns The right-hand side of the matching equation
 */
double SchroedingerSolverKronigPenney::get_rhs() const
{
    const auto L   = _l_w + _l_b; // Total length of period (well + barrier) [m]
    const auto rhs = cos(_k * L); // Right-hand side of matching equation
    return rhs;
}

/**
 * \brief Finds the left-hand side of the matching equation
 *
 * \param[in] E Energy [J]
 *
 * \returns the left-hand side of the matching equation
 */
double SchroedingerSolverKronigPenney::get_lhs(const double E) const
{
    // Wave vector inside well [QWWAD3, 2.155]
    const auto k_w = sqrt(2*_m_w*E)/hBar;

    double lhs = 0.0;

    if(E < _V0)
    {
        // If we're below the barrier, then the wavefunction decays
        // [QWWAD3, 2.167]
        const auto kappa = sqrt(2.0*_m_b*(_V0-E))/hBar;

        // Left-hand side of matching equation [QWWAD3, 2.169]
        lhs = cos(k_w*_l_w) * cosh(kappa*_l_b)
            - sin(k_w*_l_w) * sinh(kappa*_l_b) * (_m_b*_m_b*k_w*k_w - _m_w*_m_w*kappa*kappa)/(2.0*_m_w*_m_b*k_w*kappa);
    }
    else
    {
        // If we're above the barrier, then the wavefunction propagates
        // with [QWWAD3, 2.155]
        const auto k_b=sqrt(2.0*_m_b*(E-_V0))/hBar;

        // Left-hand side of matching equation [QWWAD3, 2.166]
        lhs = cos(k_w*_l_w) * cos(k_b*_l_b)
            - sin(k_w*_l_w) * sin(k_b*_l_b) * (_m_b*_m_b*k_w*k_w + _m_w*_m_w*k_b*k_b)/(2.0*_m_w*_m_b*k_w*k_b);
    }

    return lhs;
}

/**
 * \brief calculates the uncorrelated one particle wavefunction
 *
 * \param[in] E Energy [J]
 */
arma::vec SchroedingerSolverKronigPenney::get_wavefunction(const double E) const
{
    const auto nz = _z.size();
    arma::vec psi(nz); // wavefunction

#if 0
    // Temporarily disabled because this calculation doesn't work!
    // Imaginary unit
    const auto I = std::complex<double>(0,1);

    // Find local wave vector at this energy
    const auto k_w = sqrt(2.0*_m_w*E)/hBar; // wave vector in the well
    const auto k_b = sqrt(2.0*_m_b*std::complex<double>(E-_V0,0))/hBar; // Wave vector in barrier

    const auto L   = _l_w + _l_b; // Period length
#endif

    for (unsigned int i_z=0; i_z<nz; ++i_z)
    {
        // Temporarily disabled because this calculation doesn't work!
#if 0
        const auto z = _z[i_z];

        const auto A = std::complex<double>(1,0);
        const auto B = A * (sin(k_w*_l_w) + _m_b*k_w/(_m_w*k_b)*sin(k_b*_l_b)) /
            (exp(I*_k*L) * cos(k_b*_l_b) - cos(k_w*_l_w));
        const auto D = B;
        const auto C = _m_b*k_w / (_m_w*k_b) * A;
#endif

        if (_V[i_z] > 0) { // In barriers
            psi[i_z] = 1; //real(C*sin(k_b*z) + D*cos(k_b*z));
        } else { // In wells
            psi[i_z] = 1; //real(A*sin(k_w*z) + B*cos(k_w*z));
        }
    }

    return psi;
}

/**
 * \brief Test how well matched the system equations are for a given energy
 *
 * \param[in] energy Local energy
 * \param[in] solver The solver for which to evaluate parameters
 *
 * \details This gives zero for a perfect match
 */
double SchroedingerSolverKronigPenney::test_matching(double  energy,
                                                     void   *params)
{
    const auto *se = reinterpret_cast<SchroedingerSolverKronigPenney *>(params);

    const auto lhs = se->get_lhs(energy);
    const auto rhs = se->get_rhs();

    // Now evaluate the match between the left and right-hand sides
    const auto match = lhs - rhs;

    return match;
}

void SchroedingerSolverKronigPenney::calculate()
{
    gsl_function F;
    F.function  = &test_matching;
    F.params    = this;
    auto solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    const auto dx=1e-3*e; // arbitrarily small energy increment---0.1meV

    auto Elo=dx;    // first energy estimate

    for (unsigned int ist=0; ist < _nst_max; ++ist)
    {
        // Shift the lower estimate up past the last state we found
        if(ist > 0)
            Elo = _solutions[ist-1].get_energy() + dx;

        // Value for y=f(x) at bottom of search range
        const double y1 = GSL_FN_EVAL(&F,Elo);

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
        double y2; // Value of function at top of range
        double Ehi = Elo;
        do
        {
            Ehi+=dx;
            y2=GSL_FN_EVAL(&F, Ehi);
        }while((y1*y2>0)&&(Ehi<100 * _V0));

        double E; // Best estimate of energy of state
        gsl_root_fsolver_set(solver, &F, Elo, Ehi);
        int status = 0;

        // Improve the estimate of solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);

            if(status) {
                std::cerr << "GSL error in SchroedingerKronigPenney: " << std::endl
                          << "   Singularity in range (" << Elo << "," << Ehi << ")" << std::endl;
            }

            E   = gsl_root_fsolver_root(solver);
            Elo = gsl_root_fsolver_x_lower(solver);
            Ehi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(Elo, Ehi, 1e-9*e, 1e-6);
        }while(status == GSL_CONTINUE);

        // Stop if we've exceeded the cut-off energy
        if(_E_max_set && gsl_fcmp(E, _E_max, e*1e-12) == 1)
        {
            std::cerr << "Solver hit cut-off energy" << std::endl;
            break;
        }

        const auto psi = get_wavefunction(E);

        // Don't store the solution if it's below the minimum energy
        if(!(_E_min_set && gsl_fcmp(E, _E_min, e*1e-12) == -1))
        {
            _solutions.push_back(Eigenstate(E, _z, psi));
        }
    }

    gsl_root_fsolver_free(solver);
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

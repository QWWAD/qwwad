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
    SchroedingerSolver(arma::vec(nz*nper),
                       arma::vec(nz*nper),
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
            _z[iz + iper*nz] = iz*dz + iper*L;

            // Fill in barrier potential
            if(_z[iz] > l_w)
                _V[iz + iper*nz] = _V0;
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

arma::cx_mat SchroedingerSolverKronigPenney::get_matching_matrix(const double E) const
{
    const auto I = std::complex<double>(0,1);
    const auto L   = _l_w + _l_b; // Period length

    // Find local wave vector at this energy
    const auto k_w = sqrt(2.0*_m_w*E)/hBar; // wave vector in the well
    const std::complex<double> k_b = sqrt(2.0*_m_b*std::complex<double>(E-_V0,0))/hBar; // Wave vector in barrier

    arma::cx_mat M(4,4);
    M(0,0) = +exp(+I*_k *L);
    M(0,1) = +exp(+I*_k *L);
    M(0,2) = -exp(+I*k_b*L);
    M(0,3) = -exp(-I*k_b*L);
    M(1,0) = +k_w/_m_w * exp(+I*_k *L);
    M(1,1) = -k_w/_m_w * exp(+I*_k *L);
    M(1,2) = -k_b/_m_b * exp(+I*k_b*L);
    M(1,3) = +k_b/_m_b * exp(-I*k_b*L);
    M(2,0) = +exp(+I*k_w*_l_w);
    M(2,1) = +exp(-I*k_w*_l_w);
    M(2,2) = -exp(+I*k_b*_l_w);
    M(2,3) = -exp(-I*k_b*_l_w);
    M(3,0) = +k_w/_m_w * exp(+I*k_w*_l_w);
    M(3,1) = -k_w/_m_w * exp(-I*k_w*_l_w);
    M(3,2) = -k_b/_m_b * exp(+I*k_b*_l_w);
    M(3,3) = +k_b/_m_b * exp(-I*k_b*_l_w);

    return M;
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

    // Normalisation constants
    const auto I = std::complex<double>(0,1);

    // Find local wave vector at this energy
    const auto k_w = sqrt(2.0*_m_w*E)/hBar; // wave vector in the well
    const auto k_b = sqrt(2.0*_m_b*std::complex<double>(E-_V0,0))/hBar; // Wave vector in barrier

    const auto M = get_matching_matrix(E); 

    arma::cx_vec rhs(4);
    rhs(0) = 0;
    rhs(1) = 0;
    rhs(2) = 0;
    rhs(3) = 1;

    arma::cx_vec coeffs = arma::solve(M, rhs, "std");

    const auto A = coeffs(0);
    const auto B = coeffs(1);
    const auto C = coeffs(2);
    const auto D = coeffs(3);

    for (unsigned int i_z=0; i_z<nz; ++i_z)
    {
        const auto z = _z[i_z];
        const auto P = (A+B)*cos(k_w*z) + I*(A-B)*sin(k_w*z);
        const auto Q = (C+D)*cos(k_b*z) + I*(C-D)*sin(k_b*z);

        if (_V[i_z] > 0)
            psi[i_z] = 1 + abs(Q) * 0;
        else
            psi[i_z] = 1 + abs(P) * 0;
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
        double y2 = y1;
        double Ehi = Elo;
        do
        {
            Ehi+=dx;
            y2=GSL_FN_EVAL(&F, Ehi);
        }while((y1*y2>0)&&(Ehi<100 * _V0));

        auto E = (Elo + Ehi) / 2.0; // Initial estimate of energy of state
        gsl_root_fsolver_set(solver, &F, Elo, Ehi);
        int status = 0;

        // Improve the estimate of solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            E   = gsl_root_fsolver_root(solver);
            Elo = gsl_root_fsolver_x_lower(solver);
            Ehi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(Elo, Ehi, 1e-9*e, 0.000001);
        }while(status == GSL_CONTINUE);

        // Stop if we've exceeded the cut-off energy
        if(_E_cutoff_set && gsl_fcmp(E, _E_cutoff, e*1e-12) == 1)
            break;

        const auto psi = get_wavefunction(E);
        _solutions.push_back(Eigenstate(E, _z, psi));
    }

    gsl_root_fsolver_free(solver);
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

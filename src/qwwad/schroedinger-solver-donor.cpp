/**
 *  \file     schroedinger-solver-donor.cpp
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Schrodinger solver functions for donor with a hydrogenic wavefunction
 */

#include "schroedinger-solver-donor.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include "constants.h"
#include "maths-helpers.h"

namespace QWWAD
{
using namespace constants;
SchroedingerSolverDonor::SchroedingerSolverDonor(const double        m,
                                                 const decltype(_V) &V,
                                                 const decltype(_z) &z,
                                                 const double        eps,
                                                 const double        r_d,
                                                 const double        lambda,
                                                 const double        dE) :
    SchroedingerSolver(V, z, 1),
    _me(m),
    _eps(eps),
    _r_d(r_d),
    _lambda(lambda),
    _dE(dE),
    _solutions_chi()
{}

/**
 * \brief Calculates a wavefunction iteratively from left to right of structure
 *
 * \details The value of the wavefunction is taken to be zero at the point
 *          immediately to the left of the potential profile. Subsequent
 *          values are computed using QWWAD3, Eq. 5.28.
 *
 * \param[in]  E       Energy at which to compute wavefunction
 * \param[out] chi     Array to which wavefunction envelope will be written [m^{-1/2}]
 *
 * \returns The wavefunction amplitude at the point immediately to the right of the structure
 */
auto SchroedingerSolverDonor::shoot_wavefunction(const double  E,
                                                   arma::vec    &chi) const -> double
{
    const size_t nz = _z.size();
    const double dz = _z[1] - _z[0];

    chi.resize(nz);

    // boundary conditions
    chi[0] = 1;
    double chi_next = 1.0; 

    // calculate unnormalised wavefunction
    // Note that points 0 and 1 were already defined before the loop
    for(unsigned int iz = 0; iz<nz; ++iz)
    {
        // Wave function amplitude at previous point
        double chi_prev = 0.0;

        if(iz != 0)
            chi_prev = chi[iz-1];

        const double z_dash = _z[iz] - _r_d;

        const double I1=I_1(z_dash);
        const double I2=I_2(z_dash);
        const double I3=I_3(z_dash);
        const double I4=I_4(z_dash);

        const double alpha = I1;   // Coefficient of second derivative, see notes
        const double beta  = 2*I2; // Coefficient of first derivative

        // Coefficient of function
        const double gamma = I3 + _me*e*e*I4/(2.0*pi*_eps*hBar*hBar)
                             - 2.0*_me*(_V[iz]-E)*I1/(hBar*hBar);

        chi_next = ((-1.0+beta*dz/(2.0*alpha))*chi_prev
                    +(2.0-dz*dz*gamma/alpha)*chi[iz]
                   )/(1.0+beta*dz/(2.0*alpha));

        if (iz != nz - 1)
            chi[iz+1] = chi_next;
    }

    // calculate normalisation integral
    const arma::vec chi_sqr = square(chi);
    double Nchi=integral(chi_sqr,dz); // normalisation integral for chi

    /* divide unnormalised wavefunction by square root
       of normalisation integral                       */
    chi /= sqrt(Nchi);

    return chi_next/sqrt(Nchi);
}

/**
 * \brief Finds the value of the wavefunction at +infinity for a given energy.
 *
 * \details The solution to the energy occurs for chi(+infinity)=0.
 *
 * \returns The wavefunction at \f$\chi(\infty)\f$
 */
auto SchroedingerSolverDonor::chi_at_inf (double  E,
                                            void   *params) -> double
{
    const SchroedingerSolverDonor *se = reinterpret_cast<SchroedingerSolverDonor *>(params);
    arma::vec chi(se->get_z().size()); // Wavefunction envelope amplitude
    return se->shoot_wavefunction(E, chi);
}

void SchroedingerSolverDonor::calculate()
{
    _solutions_chi.clear();
    gsl_function f;
    f.function = &chi_at_inf;
    f.params   = this;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    /* initial energy estimate=minimum potential-binding energy
       of particle to free ionised dopant */
    double Elo = _V.min() - e*e/(4*pi*_eps*_lambda);

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
    double Ehi = Elo;

    double y2; // Value of f(x) at top of range
    do {
        Ehi += _dE;
        y2=GSL_FN_EVAL(&f, Ehi);
    }while(y1*y2>0);

    double E; // Best estimate of solution [J]
    gsl_root_fsolver_set(solver, &f, Elo, Ehi);
    int status = 0;

    // Improve the estimate of the solution using the Brent algorithm
    // until we hit a desired level of precision
    do
    {
        status = gsl_root_fsolver_iterate(solver);

        if(status) {
                std::cerr << "GSL error in SchroedingerSolverDonor: " << std::endl
                          << "   Singularity in range (" << Elo << "," << Ehi << ")" << std::endl;
        }

        E   = gsl_root_fsolver_root(solver);
        Elo = gsl_root_fsolver_x_lower(solver);
        Ehi = gsl_root_fsolver_x_upper(solver);
        status = gsl_root_test_interval(Elo, Ehi, 1e-12*e, 0);
    }while(status == GSL_CONTINUE);

    // Stop if we've exceeded the cut-off energy
    if(gsl_fcmp(E, _V.max()+e, e*1e-12) == 1)
        throw std::runtime_error("Energy exceeded Vmax");

    arma::vec chi(_z.size());
    const auto chi_inf = shoot_wavefunction(E, chi);
    _solutions_chi.emplace_back(E, _z, chi);

    calculate_psi_from_chi(); // Finally, compute the complete solution

    // Check that wavefunction is tightly bound
    // TODO: Implement a better check
    if(gsl_fcmp(fabs(chi_inf), 0, 1) == 1)
        throw "Warning: Wavefunction is not tightly bound";
}

/**
 * Get the solutions to the Schroedinger equation, excluding the hydrogenic component
 *
 * \details The solutions are computed on the first call to this function, but
 *          subsequent calls just recall the values and are hence much faster.
 */
auto SchroedingerSolverDonor::get_solutions_chi(const bool convert_to_meV) -> std::vector<Eigenstate>
{
    // Only calculate if we haven't done so yet
    if(_solutions_chi.empty())
    {
        calculate();

        // Delete any states that are out of the desired energy range
        // Ideally, sub-classes should never compute anything outside this
        // range!
        for(auto it = _solutions_chi.begin(); it != _solutions_chi.end(); ++it)
        {
            if ((_E_max_set && gsl_fcmp(it->get_energy(), _E_max, e*1e-12) ==  1) ||
                (_E_min_set && gsl_fcmp(it->get_energy(), _E_min, e*1e-12) == -1)) {
                _solutions_chi.erase(it);
            }
        }
    }

    if(convert_to_meV)
    {
        std::vector<Eigenstate> sol_meV;

        for(auto sol_J : _solutions_chi)
        {
            const auto E   = sol_J.get_energy();
            const auto z   = sol_J.get_position_samples();
            const auto psi = sol_J.get_wavefunction_samples();
            sol_meV.emplace_back(E, z, psi);
        }

        return sol_meV;
    }
    else
        return _solutions_chi;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

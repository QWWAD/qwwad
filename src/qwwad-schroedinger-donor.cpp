/**
 *  \file     qwwad-schroedinger-donor.cpp
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Schrodinger solver functions for donor with a hydrogenic wavefunction
 */

#include "qwwad-schroedinger-donor.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include "qclsim-constants.h"

namespace Leeds
{
using namespace constants;

SchroedingerSolverDonor::SchroedingerSolverDonor(const double                 m,
                                                 const std::valarray<double> &V,
                                                 const std::valarray<double> &z,
                                                 const double                 eps,
                                                 const double                 r_d,
                                                 const double                 lambda,
                                                 const double                 dE,
                                                 const unsigned int           nst_max) :
    SchroedingerSolver(V, z, nst_max),
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
 * \param[out] psi     Array to which complete wavefunction will be written [m^{-1/2}]
 * \param[out] chi     Array to which wavefunction envelope will be written [m^{-1/2}]
 *
 * \returns The wavefunction amplitude at the point immediately to the right of the structure
 */
double SchroedingerSolverDonor::shoot_wavefunction(const double           E,
                                                   std::valarray<double> &chi) const
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
    double Nchi=integral(pow(chi,2.0),dz); // normalisation integral for chi

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
double SchroedingerSolverDonor::chi_at_inf (double  E,
                                            void   *params)
{
    const SchroedingerSolverDonor *se = reinterpret_cast<SchroedingerSolverDonor *>(params);
    std::valarray<double> chi(se->get_z().size()); // Wavefunction envelope amplitude
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
    if(gsl_fcmp(E, _V.max(), e*1e-12) == 1)
        throw "Exceeded Vmax";

    std::valarray<double> chi(_z.size());
    const double chi_inf = shoot_wavefunction(E, chi);
    _solutions_chi.push_back(State(E,chi));

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
std::vector<State> SchroedingerSolverDonor::get_solutions_chi(const bool convert_to_meV)
{
    // Only calculate if we haven't done so yet
    if(_solutions_chi.empty())
    {
        calculate();

        // Delete any states that are out of the desired energy range
        // Ideally, sub-classes should never compute anything outside this
        // range!
        while(_E_cutoff_set && !_solutions_chi.empty() && gsl_fcmp(_solutions_chi.back().get_E(), _E_cutoff, e*1e-12) == 1)
            _solutions_chi.pop_back();

        // Normalise wavefunctions
        for(std::vector<State>::iterator ist = _solutions_chi.begin(); ist != _solutions_chi.end(); ++ist)
            ist->normalise(_z);
    }

    if(convert_to_meV)
    {
        std::vector<State> sol_meV;

        for(std::vector<State>::iterator sol_J = _solutions_chi.begin(); sol_J != _solutions_chi.end(); ++sol_J)
            sol_meV.push_back(State(sol_J->get_E()*1000/e, sol_J->psi_array()));

        return sol_meV;
    }
    else
        return _solutions_chi;
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

/** 
 * \file   fermi.cpp
 * \brief  Find the Fermi energy for a set of subbands
 * \author Alex Valavanis, University of Leeds
 * \date   2012-03-14
 * \todo   Create a program that wraps this library and outputs file
 */

#include <cstdlib>
#include <error.h>
#include <stdio.h>
#include "qclsim-constants.h"
#include "fermi.h"
#include "material_library.h"
#include "condbandpotential.h"

namespace Leeds {
using namespace constants;

/**
 * \brief Fermi occupation probability at specified kinetic energy
 *
 * \param E_F Fermi energy, relative to subband minimum [J]
 * \param Ek  Kinetic energy of electron [J]
 * \param T   temperature [K]
 *
 * \returns Fermi occupation number
 */
double f_FD(const double E_F, const double E, const double Te)
{
    return 1.0/(exp((E-E_F)/(kB*Te)) + 1.0);
}

/**
 * \brief Fermi ionisation probability with degeneracy of 2
 *
 * \param E_F Fermi energy [J]
 * \param Ed  Kinetic energy of electron [J]
 * \param T   temperature [K]
 *
 * \returns Fermi occupation number
 */
double f_FD_ionised(const double E_F, const double Ed, const double Te)
{
    return 1.0/(0.5*exp((Ed-E_F)/(kB*Te)) + 1.0);
}

/**
 * Find total population of a subband with a known Fermi energy
 * 
 * \param md  Density-of-states mass [kg]
 * \param E_F Quasi-Fermi energy, relative to subband minimum [J]
 * \param Te  Temperature of electron distribution [K]
 *
 * \returns Subband population [m^{-2}]
 */
double find_pop(const double md, const double E_F, const double Te)
{
    // Density of states in 2D system
    double rho=md/(pi*hBar*hBar);

    /* Solve Fermi integral [Harrison, 2005 (eq 2.51)] */
    /* 2011-02-23: Simplified the integral a bit, by 
     * noting that the log term can be factorised */
    double y=E_F/(kB*Te); // Just a substitution to tidy the maths
    double int_f_FD=kB*Te*gsl_log1p(exp(y)); /* Solve Fermi integral */

    return rho*int_f_FD;
}

/** 
 * Find quasi-Fermi energy for a single subband with known population and temperature
 * \param mat Material specification
 * \param N   Population density of system [m^{-2}]
 * \param Te  Temperature of carrier distribution [K]
 *
 * \returns The Fermi energy for the subband relative to the subband minimum [J]
 *
 * \todo Use dos mass calculated by Subband class.
 *       In fact, should the fermi functions all just be member functions of that
 *       class?
 */
double find_fermi(const double Esb, const double m, const double N, const double Te)
{
    const double E_F = Esb + kB*Te * log(gsl_expm1((N*pi*hBar*hBar)/(m*kB*Te)));

    return E_F;
}


#if 0
/** 
 * Find Fermi energy for an entire 2D system with many subbands, a known total population and temperature
 * \param mat Material specification
 * \param N   Population density of system [m^{-2}]
 * \param Te  Temperature of carrier distribution [K]
 * \param E   Array of subband minima [J]
 *
 * \returns The Fermi energy for the entire system [J]
 */
double find_fermi_global(const MaterialSpecification *mat,
                         const double                 N,
                         const double                 Te,
                         const std::valarray<double> &E)
{
    const size_t nst = E.size();

    // Set limits for search [J]
    double E_min=E[0]-100.0*kB*Te;
    double E_max=E[nst-1]+100.0*kB*Te;

    // Find bisector of the limits [J]
    double E_mid=(E_min+E_max)/2.0;

    // TODO: Should calculate band-edge parameters for a given substrate
    //       and temperature.
    const BandEdge band_edge(*mat, *mat);
    const double _md = band_edge.get_md(); // Density-of-states mass [kg]

    // Find the signs of ∫f(E_min) - N at the endpoints
    const int sign_min=GSL_SIGN(find_pop(_md, E_min, Te) -  N);
    const int sign_max=GSL_SIGN(find_pop(_md, E_max, Te) -  N);

    // We can solve the Fermi integral only if there is a sign-change 
    // between the limits
    const bool solvable=!(sign_min==sign_max);

    if(solvable)
    {
        // Solve iteratively using linear-bisection to a precision of
        // 0.01 micro-eV
        // TODO: Make precision configurable
        while(fabs(E_max-E_min) > 1e-8*e)
        {
            double total_population = 0.0;

            for(unsigned int ist = 0; ist < nst; ist++)
                total_population += find_pop(_md, E_mid - E[ist], Te);

            const int sign_mid=GSL_SIGN(total_population - N);

            if(sign_mid==sign_min) E_min=E_mid;
            else E_max=E_mid;

            E_mid=(E_max+E_min)/2.0;
        }
    }
    else error(EXIT_FAILURE, 0, "ERROR: No quasi-Fermi energy in range.");

    return E_mid;
}
#endif
} // namespace Leeds

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

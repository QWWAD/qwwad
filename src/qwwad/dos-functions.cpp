/**
 * \file dos-functions.cpp Functions for calculating density of states
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "dos-functions.h"
#include "constants.h"
#include <cmath>

namespace QWWAD {
using namespace constants;

/**
 * Calculate bulk density of states
 *
 * \param[in] mass   The effective mass of the carrier at the band edge [kg]
 * \param[in] energy The energy of the carrier above the band edge [J]
 * \param[in] V      Band edge potential [J]
 * \param[in] alpha  The nonparabolicity factor [1/J]
 *
 * \returns The density of states [J^{-1}m^{-3}]
 *
 * \details The effective mass is assumed to take the form
 *          \f[
 *          m^*(E) = m_0[1 + \alpha(E-V)]
 *          \f]
 */
auto calculate_dos_3D(const double mass,
                        const double energy,
                        const double V,
                        const double alpha) noexcept -> double
{
    double rho = 0.0;

    // Only bother calculating if we're above the band edge
    if(gsl_fcmp(energy,V,energy/1e12)==1)
    {
        const double E = energy - V; // Express energy relative to band edge

        // Density-of-states mass for an excited state in bulk
        const auto m_d = mass * std::cbrt((1.0+alpha*E) * gsl_pow_2(1.0+2*alpha*E));

        // Bulk density of states [QWWAD3, Eq. 2.40]
        rho = 1/(2.0*gsl_pow_2(pi)) * gsl_pow_3(sqrt(2.0*m_d/gsl_pow_2(hBar))) * sqrt(E);
    }

    return rho;
}

/**
 * Calculate 2D density of states for a quantum well
 *
 * \param[in] mass       The effective mass of the carrier [kg]
 * \param[in] E_carrier  The energy of the carrier [J]
 * \param[in] E_subbands The energies of the minima of each subband in the well [J]
 * \param[in] V          Band edge potential [J]
 * \param[in] alpha      The nonparabolicity factor [1/J]
 *
 * \returns The density of states [J^{-1}m^{-2}]
 *
 * \details The effective mass is assumed to take the form
 *          \f[
 *          m^*(E) = m_0[1 + \alpha(E-V)]
 *          \f]
 */
auto calculate_dos_2D(const double                 mass,
                        const double                 E_carrier,
                        const std::valarray<double> &E_subbands,
                        const double                 V,
                        const double                 alpha) noexcept -> double
{
    // Density of states in a single subband [QWWAD3, Eq. 2.46]
    const double m_d = mass * (1 + 2*alpha*(E_carrier-V));
    const double dos_1sb = m_d/(pi*gsl_pow_2(hBar));

    double dos_total = 0; // Total dos over all occupied subbands
    const size_t nsb = E_subbands.size();

    // Loop over subbands
    for(unsigned int isb = 0; isb < nsb; ++isb) {
        // Increment dos whenever we enter a new subband
        // [QWWAD3, Eq. 2.47]
        if(E_carrier > E_subbands[isb]) {
            dos_total += dos_1sb;
        }
    }

    return dos_total;
}

/**
 * Calculate 1D density of states for a quantum wire
 *
 * \param[in] mass       The effective mass of the carrier [kg]
 * \param[in] E_carrier  The energy of the carrier [J]
 * \param[in] E_subbands The energies of the minima of each subband in the well [J]
 *
 * \returns The density of states [J^{-1}m^{-1}]
 *
 * \details It is assumed that carriers have a parabolic dispersion
 */
auto calculate_dos_1D(const double                 mass,
                        const double                 E_carrier,
                        const std::valarray<double> &E_subbands) noexcept -> double
{
    double dos_total = 0; // Total dos over all occupied subbands
    const size_t nsb = E_subbands.size();

    // Loop over subbands
    for(unsigned int isb = 0; isb < nsb; ++isb) {
        // Increment dos whenever we enter a new subband
        if(E_carrier > E_subbands[isb]) {
            dos_total += sqrt(2*mass)/hBar/(pi*sqrt(E_carrier - E_subbands[isb]));
        }
    }

    return dos_total;
}
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

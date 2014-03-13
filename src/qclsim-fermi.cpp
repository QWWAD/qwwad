/** 
 * \file   qclsim-fermi.cpp
 * \brief  Find the Fermi energy for a set of subbands
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "qclsim-constants.h"
#include "qclsim-fermi.h"
#include <stdexcept>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_fermi_dirac.h>

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
 * \brief Find total population of a subband with a known Fermi energy
 * 
 * \param Esb Energy of the subband minimum [J]
 * \param E_F Quasi-Fermi energy on the same absolute scale as the Fermi energy [J]
 * \param m0  Band-edge effective mass [kg]
 * \param Te  Temperature of electron distribution [K]
 *
 * \returns Subband population [m^{-2}]
 */
double find_pop(const double Esb,
                const double E_F,
                const double m0,
                const double Te,
                const double alpha,
                const double V)
{
    double N = 0; // Population to output

    // Density of states in 2D system with parabolic dispersion
    const double rho_p = m0/(pi*hBar*hBar);

    const double x = (E_F - Esb)/(kB*Te); // Substitution to simplify the expression

    // Use the parabolic (simple) solution if possible
    if(gsl_fcmp(alpha,0,1e-6) == 0)
    {
        // Solve Fermi integral (eq 2.66, QWWAD4)
        N = rho_p*kB*Te*gsl_sf_fermi_dirac_0(x);
    }
    else
    {
        // Full non-parabolic solution (eq 2.69, QWWAD4)
        N = rho_p * kB*Te *(
                (1.0 + 2.0 * alpha * (Esb-V)) * gsl_sf_fermi_dirac_0(x)
                + 2*alpha*kB*Te * gsl_sf_fermi_dirac_1(x)
                );
    }

    return N;
}

/** 
 * \brief Find quasi-Fermi energy for a single subband with known population and temperature
 *
 * \param[in] Esb Energy of the subband minimum [J]
 * \param[in] m   Mass of carriers [kg]
 * \param[in] N   Population density of system [m^{-2}]
 * \param[in] Te  Temperature of carrier distribution [K]
 *
 * \returns The Fermi energy for the subband [J]
 *
 * \details The Fermi energy is calculated using equation 2.58, QWWAD4
 *          \f[
 *            E_{\text{F}} = E_{\text{min}} + kT\ln{\left[\mbox{e}^{\left(\frac{N\pi\hbar^2}{m^*kT}\right)}-1\right]}
 *          \f]
 *          which assumes that the carriers can spread to any energy above the subband
 *          minimum.
 *
 * \todo Use dos mass calculated by Subband class.
 *       In fact, should the fermi functions all just be member functions of that
 *       class?
 */
double find_fermi(const double Esb,
                  const double m,
                  const double N,
                  const double Te,
                  const double alpha,
                  const double V)
{
    double E_F = 0;

    // Use the analytical form if possible
    if(gsl_fcmp(alpha, 0.0, 1.0e-6) == 0)
    {
        // Eq. 2.75, QWWAD4
        E_F = Esb + kB*Te * log(gsl_expm1((N*pi*hBar*hBar)/(m*kB*Te)));
    }
    else
    {
        // Set limits for search [J]
        double E_min=Esb-100.0*kB*Te;
        double E_max=Esb+100.0*kB*Te;

        // Find bisector of the limits [J]
        double E_mid=(E_min+E_max)/2.0;

        // Find the signs of ∫f(E_min) - N at the endpoints
        const int sign_min=GSL_SIGN(find_pop(Esb, E_min, m, Te, alpha, V) -  N);
        const int sign_max=GSL_SIGN(find_pop(Esb, E_max, m, Te, alpha, V) -  N);

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
                double pop = find_pop(Esb, E_mid, m, Te, alpha, V);

                const int sign_mid=GSL_SIGN(pop - N);

                if(sign_mid==sign_min) E_min=E_mid;
                else E_max=E_mid;

                E_mid=(E_max+E_min)/2.0;
            }
        }
        else throw std::runtime_error("No quasi-Fermi energy in range.");

        E_F = E_mid;
    }

    return E_F;
}

/** 
 * \brief Find Fermi energy for an entire 2D system with many subbands, a known total population and temperature
 *
 * \param m   Mass of carrier [kg]
 * \param N   Population density of system [m^{-2}]
 * \param Te  Temperature of carrier distribution [K]
 * \param E   Array of subband minima [J]
 *
 * \returns The Fermi energy for the entire system [J]
 */
double find_fermi_global(const double                 m,
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

    // Find the signs of ∫f(E_min) - N at the endpoints
    const int sign_min=GSL_SIGN(find_pop(E[0], E_min, m, Te) -  N);
    const int sign_max=GSL_SIGN(find_pop(E[0], E_max, m, Te) -  N);

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
                total_population += find_pop(E[ist], E_mid, m, Te);

            const int sign_mid=GSL_SIGN(total_population - N);

            if(sign_mid==sign_min) E_min=E_mid;
            else E_max=E_mid;

            E_mid=(E_max+E_min)/2.0;
        }
    }
    else throw std::runtime_error("No quasi-Fermi energy in range.");

    return E_mid;
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

/** 
 * \file   qclsim-fermi.cpp
 * \brief  Find the Fermi energy for a set of subbands
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "qwwad/constants.h"
#include "qclsim-fermi.h"
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_fermi_dirac.h>

namespace QWWAD
{
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
 * \param V   Energy of the band edge [J]
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

    // In case of underflow, return a tiny value
    if(gsl_fcmp(x,-700,1e-6) == -1)
            return 1;

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
 * \brief Parameters used for checking the population
 */
struct pop_params {
    double Esb;     ///< Subband minimum [J]
    double m0;      ///< Band-edge effective mass [kg]
    double Te;      ///< Temperature of carrier distribution [K]
    double alpha;   ///< Nonparabolicity [1/J]
    double V;       ///< Band-edge [J]
    double N;       ///< Population we expect to find [m^{-3}]
};

/**
 * \brief Find the error in population for a given Fermi energy
 */
double find_pop_error(double E_F, void *params)
{
    const pop_params *p = reinterpret_cast<pop_params *>(params);
    return find_pop(p->Esb, E_F, p->m0, p->Te, p->alpha, p->V) - p->N;
}

/**
 * \brief Parameters used for checking the population of multiple subbands
 */
struct total_pop_params {
    std::valarray<double> Esb; ///< Subband minimum [J]
    double m0;                 ///< Band-edge effective mass [kg]
    double Te;                 ///< Temperature of carrier distribution [K]
    double alpha;              ///< Nonparabolicity [1/J]
    double V;                  ///< Band-edge [J]
    double N;                  ///< Population we expect to find [m^{-3}]
};

/**
 * \brief Find the error in the total population for multiple subbands
 */
double find_total_pop_error(double E_F, void *params)
{
    const total_pop_params *p = reinterpret_cast<total_pop_params *>(params);

    // Find total population in each subband, using the same global Fermi
    // energy
    double N_total = 0.0;
    const size_t nst = p->Esb.size();

    for(unsigned int ist = 0; ist < nst; ist++)
        N_total += find_pop(p->Esb[ist], E_F, p->m0, p->Te, p->alpha, p->V);

    return N_total - p->N;
}


/** 
 * \brief Find quasi-Fermi energy for a single subband with known population and temperature
 *
 * \param[in] Esb   Energy of the subband minimum [J]
 * \param[in] m     Mass of carriers [kg]
 * \param[in] N     Population density of system [m^{-2}]
 * \param[in] Te    Temperature of carrier distribution [K]
 * \param[in] alpha Nonparabolicity parameter [1/J]
 * \param[in] V     Band-edge [J]
 *
 * \returns The Fermi energy for the subband [J]
 *
 * \details For parabolic bands, the Fermi energy is calculated using equation 2.58, QWWAD4
 *          \f[
 *            E_{\text{F}} = E_{\text{min}} + kT\ln{\left[\mbox{e}^{\left(\frac{N\pi\hbar^2}{m^*kT}\right)}-1\right]}
 *          \f]
 *          which assumes that the carriers can spread to any energy above the subband
 *          minimum.
 *
 *          For nonparabolic bands, the Fermi energy is found numerically
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

        // Set parameters for search
        pop_params params = {Esb, m, Te, alpha, V, N};

        // Find the signs of ∫f(E_min) - N at the endpoints
        const int sign_min=GSL_SIGN(find_pop_error(E_min, &params));
        const int sign_max=GSL_SIGN(find_pop_error(E_max, &params));

        // We can solve the Fermi integral only if there is a sign-change 
        // between the limits
        const bool solvable=!(sign_min==sign_max);

        if(solvable)
        {
            gsl_function F;
            F.function = &find_pop_error;
            F.params   = &params;
            gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
            gsl_root_fsolver_set(solver, &F, E_min, E_max);

            int status;
            // Solve iteratively
            do
            {
                status = gsl_root_fsolver_iterate(solver);
                E_F   = gsl_root_fsolver_root(solver);
                E_min = gsl_root_fsolver_x_lower(solver);
                E_max = gsl_root_fsolver_x_upper(solver);
                status = gsl_root_test_interval(E_min, E_max, 1e-8*e, 0);
            }while(status == GSL_CONTINUE);

            gsl_root_fsolver_free(solver);
        }
        else throw std::runtime_error("No quasi-Fermi energy in range.");
    }

    return E_F;
}

/** 
 * \brief Find Fermi energy for an entire 2D system with many subbands, a known total population and temperature
 *
 * \param states Array of states in system
 * \param m0     Mass of carriers at band edge [kg]
 * \param N      Population density of system [m^{-2}]
 * \param Te     Temperature of carrier distribution [K]
 * \param alpha  Nonparabolicity parameter [1/J]
 * \param V      Band-edge [J]
 *
 * \returns The Fermi energy for the entire system [J]
 */
double find_fermi_global(const std::vector<Eigenstate> &states,
                         const double                   m0,
                         const double                   N,
                         const double                   Te,
                         const double                   alpha,
                         const double                   V)
{
    std::valarray<double> E(states.size());

    for(unsigned int ist = 0; ist < E.size(); ++ist)
        E[ist] = states[ist].get_energy();

    return find_fermi_global(E, m0, N, Te, alpha, V);
}

/** 
 * \brief Find Fermi energy for an entire 2D system with many subbands, a known total population and temperature
 *
 * \param Esb   Array of subband minima [J]
 * \param m0    Mass of carriers at band edge [kg]
 * \param N     Population density of system [m^{-2}]
 * \param Te    Temperature of carrier distribution [K]
 * \param alpha Nonparabolicity parameter [1/J]
 * \param V     Band-edge [J]
 *
 * \returns The Fermi energy for the entire system [J]
 */
double find_fermi_global(const std::valarray<double> &Esb,
                         const double                 m0,
                         const double                 N,
                         const double                 Te,
                         const double                 alpha,
                         const double                 V)
{
    const size_t nst = Esb.size();

    // Set limits for search [J]
    double E_min=Esb[0]-100.0*kB*Te;
    double E_max=Esb[nst-1]+500.0*kB*Te;

    // Find bisector of the limits [J]
    double E_F=(E_min+E_max)/2.0;

    // Set parameters for search
    total_pop_params params = {Esb, m0, Te, alpha, V, N};

    // Find the signs of ∫f(E_min) - N at the endpoints
    const int sign_min=GSL_SIGN(find_total_pop_error(E_min, &params));
    const int sign_max=GSL_SIGN(find_total_pop_error(E_max, &params));

    // We can solve the Fermi integral only if there is a sign-change 
    // between the limits
    const bool solvable=!(sign_min==sign_max);

    if(solvable)
    {
        gsl_function F;
        F.function = &find_total_pop_error;
        F.params   = &params;
        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        gsl_root_fsolver_set(solver, &F, E_min, E_max);

        int status;
        // Solve iteratively
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            E_F   = gsl_root_fsolver_root(solver);
            E_min = gsl_root_fsolver_x_lower(solver);
            E_max = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(E_min, E_max, 1e-8*e, 0);
        }while(status == GSL_CONTINUE);

        gsl_root_fsolver_free(solver);
    }
    else throw std::runtime_error("No quasi-Fermi energy in range.");

    return E_F;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

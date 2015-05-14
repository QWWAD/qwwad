/**
 * \file   hup.cpp
 * \brief  Heisenberg's Uncertainty Principle
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program calculates the uncertainty in the position and the
 *          momentum of any user supplied wave function.  The data is written 
 *          to the standard output.
 */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qclsim-linalg.h"
#include "wf_options.h"
#include "qwwad/constants.h"

using namespace Leeds;
using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
WfOptions configure_options(int argc, char* argv[])
{
    WfOptions opt(WF_OPTION_MODE_IN);
    
    std::string summary("Compute the uncertainty relation for a given state.");

    opt.add_size_option("state",  1, "Number of state to analyse.");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc, char *argv[])
{
    WfOptions opt = configure_options(argc, argv);

    unsigned int state = opt.get_size_option("state"); // Principal quantum number of state to analyze

    std::vector<State> all_states = State::read_from_file(opt.get_energy_input_path(),
                                                          opt.get_wf_input_prefix(),
                                                          opt.get_wf_input_ext(),
                                                          1000.0/e,
                                                          true);

    // TODO: read z-coordinates in from the ground state
    // There should probably be an option in the State class for this
    std::valarray<double> z;
    std::valarray<double> _psi;
    read_table("wf_e1.r", z, _psi);
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    std::valarray<double> psi = all_states.at(state-1).psi_array();

    std::valarray<double> d_psi_dz(0.0, nz);   // Derivative of wavefunction
    std::valarray<double> d2_psi_dz2(0.0, nz); // 2nd Derivative of wavefunction

    // Note that we can take the end points as zero, since this is
    // guaranteed for any valid wavefunction
    for(unsigned int i=1;i<nz-1;i++)
    {
        d_psi_dz[i]   = (psi[i+1] - psi[i-1])/(2*dz);
        d2_psi_dz2[i] = (psi[i+1] - 2*psi[i] + psi[i-1])/(dz*dz);
    }

    const std::valarray<double> ev_z_integrand_dz = psi*z*psi;
    const double ev_z    = integral(ev_z_integrand_dz, dz);       // Expectation position [m]

    const std::valarray<double> ev_zsqr_integrand_dz = psi*z*z*psi;
    const double ev_zsqr = integral(ev_zsqr_integrand_dz, dz);     // Expectation for z*z [m^2]

    const std::valarray<double> ev_p_integrand_dz = -psi*d_psi_dz;
    const double ev_p    = integral(ev_p_integrand_dz, dz);   // Expectation momentum [relative to i hBar]

    const std::valarray<double> ev_psqr_integrand_dz = -psi*d2_psi_dz2;
    const double ev_psqr = integral(ev_psqr_integrand_dz, dz); // Expectation for p*p [relative to hBar^2]

    // Find uncertainty in position and momentum
    const double Delta_z=sqrt(ev_zsqr-gsl_pow_2(ev_z));
    const double Delta_p=sqrt(ev_psqr-gsl_pow_2(ev_p));

    printf("<z>\t\t\t\t\t%20.17le\n",ev_z);
    printf("<z^2>\t\t\t\t\t%20.17le\n",ev_zsqr);
    printf("Delta_z=sqrt(<z^2>-<z>^2)\t\t%20.17le\n",Delta_z);

    printf("\n");

    printf("<p>/hbar\t\t\t\t%20.17le i\n",ev_p);
    printf("<p^2>/sqr(hbar)\t\t\t\t%20.17le\n",ev_psqr);
    printf("Delta_p/hbar=sqrt(<p^2>-<p>^2)/hbar\t%20.17le\n",Delta_p);

    printf("\n");

    printf("Delta_z*Delta_p\t\t\t\t%20.17le hbar\n",Delta_z*Delta_p);

    return EXIT_SUCCESS;
}

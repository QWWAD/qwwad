/*====================================================================
                  hup Heisenberg's Uncertainty Principle
  ====================================================================*/

/* This program calculates the uncertainty in the position and the
   momentum of any user supplied wave function.  The data is written 
   to the standard output.

   Input files:       wf_ps.r		p=particle (e, h or l)
					s=state

   Output files:    

   Paul Harrison, April 1998						*/
                                                               

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "qclsim-linalg.h"
#include "wf_options.h"
#include "qclsim-constants.h"

using namespace Leeds;
using namespace Leeds::constants;

/**
 * Configure command-line options for the program
 */
WfOptions configure_options(int argc, char* argv[])
{
    WfOptions opt;
    opt.add_size_option("state",  1, "Number of state to analyse.");

    std::string summary("Compute the uncertainty relation for a given state.");
    std::string details("Default input files:\n"
                        "  'E*.r'    \tEnergy of each state:\n"
                        "            \tCOLUMN 1: state index.\n"
                        "            \tCOLUMN 2: energy [meV].\n"
                        "  'wf_*i.r' \tWave function amplitude at each position\n"
                        "            \tCOLUMN 1: position [m].\n"
                        "            \tCOLUMN 2: wave function amplitude [m^{-1/2}].\n"
                        "\n"
                        "\tIn each case, the '*' is replaced by the particle ID and the 'i' is replaced by the number of the state.\n"
                        "\n"
                        "Default output files:\n"
                        "\n"
                        "Examples:\n"
                        "   Calculate uncertainty in <z><p> for state 3:\n\n"
                        "   hup --state 3\n");

    opt.add_prog_specific_options_and_parse(argc, argv, summary, details);

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
read_table_xy("wf_e1.r", z, _psi);
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

const double ev_z    = trapz(psi*z*psi, dz);       // Expectation position [m]
const double ev_zsqr = trapz(psi*z*z*psi, dz);     // Expectation for z*z [m^2]
const double ev_p    = trapz(-psi*d_psi_dz, dz);   // Expectation momentum [relative to i hBar]
const double ev_psqr = trapz(-psi*d2_psi_dz2, dz); // Expectation for p*p [relative to hBar^2]

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

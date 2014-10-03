/**
 * \file   tsb.cpp
 * \brief  Calculate transmission coefficient for single barrier
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Find the transmission coefficient through a single tunnelling barrier.");

    opt.add_numeric_option("barrier-width,L",    100,   "Width of barrier [angstrom].");
    opt.add_numeric_option("mass,m",             0.067, "Effective mass in barrier (relative to that of a free electron). "
                                                        "This is assumed to be constant throughout the whole system.");
    opt.add_numeric_option("potential",           100,  "Barrier potential [meV]");
    opt.add_numeric_option("energy-step,d",        0.1, "Energy step [meV]");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double dE = opt.get_numeric_option("energy-step") * 1e-3*e;   // [J]
    const double L  = opt.get_numeric_option("barrier-width") * 1e-10;
    const double m  = opt.get_numeric_option("mass") * me;
    const double V  = opt.get_numeric_option("potential") * e / 1000;

    const double E_cutoff = V * 10; // Cut-off energy for plot
    const size_t nE = floor(E_cutoff/dE); // Number of points in plot

    std::valarray<double> E(nE); // Array of energies
    std::valarray<double> T(nE); // Array of transmission coefficients

    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE; // Find energy
        const double k=sqrt(2*m*E[iE])/hBar; // Wave-vector in `well' material

        if(gsl_fcmp(E[iE], V, dE/1000) == -1) // If E < V
        {
            const double K = sqrt(2*m*(V-E[iE]))/hBar; // Decay constant in barrier
            T[iE] = 1/(1+gsl_pow_2((k*k+K*K)/(2*k*K) * sinh(K*L))); // [QWWAD4, 2.199]
        }
        else // if E > V
        {
            const double kdash = sqrt(2*m*(E[iE]-V))/hBar; // Wave-vector above barrier
            T[iE] = 1/(1+gsl_pow_2((k*k-kdash*kdash)/(2*k*kdash) * sin(kdash*L))); // [QWWAD4, 2.201]
        }
    }

    // Rescale to meV for output
    E/=(1e-3*e);
    write_table_xy("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

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
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string summary("Find the transmission coefficient through a single tunnelling barrier.");

    opt.add_option<double>("barrierwidth,L",    100,  "Width of barrier [angstrom].");
    opt.add_option<double>("mass,m",           0.067, "Effective mass (relative to that of a free electron). "
                                                      "This is assumed to be constant throughout the whole system.");
    opt.add_option<double>("barrierpotential",  100,  "Barrier potential [meV]");
    opt.add_option<double>("dE,d",              0.1,  "Energy step [meV]");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto dE = opt.get_option<double>("dE") * 1e-3*e;   // [J]
    const auto L  = opt.get_option<double>("barrierwidth") * 1e-10;
    const auto m  = opt.get_option<double>("mass") * me;
    const auto V  = opt.get_option<double>("barrierpotential") * e / 1000;

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
            T[iE] = 1.0/(1.0 + gsl_pow_2((k*k+K*K)/(2*k*K) * sinh(K*L))); // [QWWAD4, 2.199]
        }
        else if(gsl_fcmp(E[iE], V, dE/1000) == 0) // If E = V
        {
            // Use analytical solution, which is obtained from L'Hopital's rule
            T[iE] = 4.0/(4.0 + k*k*L*L);
        }
        else // if E > V
        {
            const double kdash = sqrt(2*m*(E[iE]-V))/hBar; // Wave-vector above barrier
            T[iE] = 1/(1+gsl_pow_2((k*k-kdash*kdash)/(2*k*kdash) * sin(kdash*L))); // [QWWAD4, 2.201]
        }
    }

    // Rescale to meV for output
    E/=(1e-3*e);

    try {
        write_table("T.r", E, T);
    } catch (std::runtime_error &e) {
        std::cerr << "Error writing to file" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

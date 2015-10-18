/**
 * \file   qwwad_tx_double_barrier.cpp
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details calculates the transmission coefficient T(E)
 *          for a double barrier structure as a function of the 
 *          incident energy E of a particle.  The barriers can be of 
 *          unequal width and the particle can have a different mass in
 *          the barrier material to the `well', the barrier heights are
 *          however equal.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <valarray>
#include "qwwad/double-barrier.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Find the transmission coefficient through a double tunnelling barrier.");

    opt.add_option<double>("leftbarrierwidth,l",      100, "Width of left barrier [angstrom].");
    opt.add_option<double>("wellwidth,b",             100, "Width of well [angstrom].");
    opt.add_option<double>("rightbarrierwidth,r",     100, "Width of right barrier [angstrom].");
    opt.add_option<double>("wellmass,m",            0.067, "Effective mass in well (relative to that of a free electron).");
    opt.add_option<double>("barriermass,n",         0.067, "Effective mass in barrier (relative to that of a free electron).");
    opt.add_option<double>("barrierpotential",        100, "Barrier potential [meV]");
    opt.add_option<double>("dE,d",                   0.01, "Energy step [meV]");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto dE  = opt.get_option<double>("dE") * 1e-3 * e;               // [J]
    const auto L1  = opt.get_option<double>("leftbarrierwidth") * 1e-10;    // [m]
    const auto L2  = opt.get_option<double>("wellwidth") * 1e-10;           // [m]
    const auto L3  = opt.get_option<double>("rightbarrierwidth") * 1e-10;   // [m]
    const auto m_w = opt.get_option<double>("wellmass") * me;               // [kg]
    const auto m_b = opt.get_option<double>("barriermass") * me;            // [kg]
    const auto V   = opt.get_option<double>("barrierpotential") * e / 1000; // [J]

    const size_t nE = floor(V/dE); // Number of points in plot

    std::valarray<double> E(nE); // Array of energies
    std::valarray<double> T(nE); // Array of transmission coefficients

    // Loop over energy
    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE;
        T[iE] = get_transmission_coefficient(E[iE], m_w, m_b, V, L1, L2, L3);
    }

    // Rescale to meV for output
    E/=(1e-3*e);
    write_table("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

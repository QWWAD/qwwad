/**
 * \file   tdb.cpp
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
#include "double-barrier.h"
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
    opt.add_numeric_option("left-barrier-width,l",    100, "Width of left barrier [angstrom].");
    opt.add_numeric_option("well-width,b",            100, "Width of well [angstrom].");
    opt.add_numeric_option("right-barrier-width,r",   100, "Width of right barrier [angstrom].");
    opt.add_numeric_option("well-mass,m",           0.067, "Effective mass in well (relative to that of a free electron).");
    opt.add_numeric_option("barrier-mass,n",        0.067, "Effective mass in barrier (relative to that of a free electron).");
    opt.add_numeric_option("potential",               100, "Barrier potential [meV]");
    opt.add_numeric_option("energy-step,d",          0.01, "Energy step [meV]");

    std::string summary("Find the transmission coefficient through a double tunnelling barrier.");
    std::string details("The following output text file is created:\n"
                        "  'T.r'    \tTransmission coefficient as a function of energy:\n"
                        "           \tCOLUMN 1: energy of incident carrier [meV].\n"
                        "           \tCOLUMN 2: transmission coefficient."
                        "\n"
                        "Examples:\n"
                        "   Compute transmission coefficient through a pair of 100-meV, 200-angstrom barriers separated by a 50-angstrom well:\n\n"
                        "   tdb --left-barrier-width 200 --right-barrier-width 200 --potential 100 --well-width 50\n"
                        "\n"
                        "   Compute transmission coefficient through a pair of 1-eV, 100-angstrom barrier, with effective mass = 0.1 m0, using resolution of 1 meV:\n\n"
                        "   tdb --potential 1000 --left-barrier-width 100 --right-barrier-width 100 --barrier-mass 0.1 --energy-step 1\n");

    opt.add_prog_specific_options_and_parse(argc, argv, summary, details);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double dE  = opt.get_numeric_option("energy-step") * 1e-3 * e;      // [J]
    const double L1  = opt.get_numeric_option("left-barrier-width") * 1e-10;  // [m]
    const double L2  = opt.get_numeric_option("well-width") * 1e-10;          // [m]
    const double L3  = opt.get_numeric_option("right-barrier-width") * 1e-10; // [m]
    const double m_w = opt.get_numeric_option("well-mass") * me;              // [kg]
    const double m_b = opt.get_numeric_option("barrier-mass") * me;           // [kg]
    const double V   = opt.get_numeric_option("potential") * e / 1000;        // [J]

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
    write_table_xy("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

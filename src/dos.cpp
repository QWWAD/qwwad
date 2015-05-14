/**
 * \file   dos.cpp
 * \brief  Density of states calculator
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program calculates the density of states for bulk (3D),
 *          quantum wells (2D) and quantum wires (1D), for a series of subband
 *          minima which are read in from the external file `Ee.r', or `Eh.r'
 */

#include <iostream>
#include "dos-functions.h"
#include "qwwad/constants.h"
#include "qwwad/file-io-deprecated.h"
#include "qwwad-options.h"
#include <cstdlib>
#include <valarray>
#include <fstream>

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Find density of states for bulk (3D), quantum wells (2D) and quantum wires (1D).");

    opt.add_numeric_option("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_numeric_option("vcb",        0.00,  "Band-edge potential [eV]");
    opt.add_numeric_option("alpha",      0.00,  "Non-parabolicity parameter [eV^{-1}]");
    opt.add_char_option   ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, heavy holes or light holes respectively.");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
};

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto p     = opt.get_char_option("particle");         // particle ID (e, h or l)
    const auto m     = opt.get_numeric_option("mass") * me;     // effective mass [kg]
    const auto V     = opt.get_numeric_option("vcb") * e;       // band_edge potential [J]
    const auto alpha = opt.get_numeric_option("alpha") / e;     // Non-parabolicity [1/J]

    const size_t n=1000; // Number of output energies

    const auto E = read_E(p); // read in subband minima [J]

    std::valarray<double> energy(n+1);   // Energies at which dos is calculated [J]
    std::valarray<double> dos_bulk(n+1); // bulk (3D) dos [J^{-1}m^{-3}]
    std::valarray<double> dos_2D(n+1);   // quantum well (2D) dos [J^{-1}m^{-2}]
    std::valarray<double> dos_1D(n+1);   // quantum wire (1D) dos [J^{-1}m^{-1}]

    // TODO: Use QWWAD Fileio instead
    std::ofstream Frho("rho.r");

    for(unsigned int ie=0;ie<=n;ie++)
    {
        energy[ie] = ie*1e-3*e; // convert meV-> J

        dos_bulk[ie] = calculate_dos_3D(m, energy[ie], V, alpha);
        dos_2D[ie]   = calculate_dos_2D(m, energy[ie], E, V, alpha);
        dos_1D[ie]   = calculate_dos_1D(m, energy[ie], E);

        Frho << energy[ie]/(1e-3*e) << " " << dos_bulk[ie] << " " << dos_2D[ie] << " " << dos_1D[ie] << std::endl;
    }

    Frho.close();
    return EXIT_SUCCESS;
} /* end main */
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :

/**
 * \file   dos.c Density of states calculator
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program calculates the density of states for bulk (3D),
 *          quantum wells (2D) and quantum wires (1D), for a series of subband
 *          minima which are read in from the external file `Ee.r', or `Eh.r'
 */

#include <iostream>
#include "dos-functions.h"
#include "qclsim-constants.h"
#include "qwwad-fileio.h"
#include "qwwad-options.h"
#include <cstdlib>
#include <valarray>
#include <fstream>

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class DOSOptions : public Options
{
    public:
        DOSOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass (relative to that of a free electron)")

                    ("band-edge", po::value<double>()->default_value(0),
                     "Band edge potential [eV]")

                    ("alpha", po::value<double>()->default_value(0),
                     "Nonparabolicity parameter [eV^{-1}]")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")
                    ;

                std::string doc("Find the density of states for bulk (3D), "
                                "quantum wells (2D) and quantum wires (1D), "
                                "for a series of subband minima which are read "
                                "in from the external file `Ee.r', or `Eh.r'."
                                "Data is written to the file `rho.r' in the format: \n"
                                "Column 1 = energy (meV)\n"
                                "Column 2 = 3D density of states (J^{-1}m^{-3})\n"
                                "Column 3 = 2D density of states (J^{-1}m^{-2})\n"
                                "Column 4 = 1D density of states (J^{-1}m^{-1})");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the effective mass [kg]
        double get_mass() const {return vm["mass"].as<double>()*me;}

        /// \returns the band edge [J]
        double get_band_edge() const {return vm["band-edge"].as<double>()*e;}

        /// \returns the nonparabolicity parameter [J^{-1}]
        double get_alpha() const {return vm["alpha"].as<double>()/e;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}
};

int main(int argc,char *argv[])
{
    DOSOptions opt(argc, argv);
    const char   p = opt.get_particle();  // particle (e, h or l)
    const double m = opt.get_mass();      // effective mass [kg]
    const double V = opt.get_band_edge(); // band-edge potential [J]
    const double alpha = opt.get_alpha(); // Nonparabolicity parameter [1/J]

    const size_t n=1000; // Number of output energies

    std::valarray<double> E = read_E(p); // read in subband minima [J]

    std::valarray<double> energy(n+1);   // Energies at which dos is calculated [J]
    std::valarray<double> dos_bulk(n+1); // bulk (3D) dos [J^{-1}m^{-3}]
    std::valarray<double> dos_2D(n+1);   // quantum well (2D) dos [J^{-1}m^{-2}]
    std::valarray<double> dos_1D(n+1);   // quantum wire (1D) dos [J^{-1}m^{-1}]

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

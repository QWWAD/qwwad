/**
 * \file   efiw.cpp Calculates eigenstates of an infinite square well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "qwwad-schroedinger-infinite-well.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    opt.add_numeric_option("width,L",    100,   "Width of quantum well [angstrom].");
    opt.add_numeric_option("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_size_option   ("nz,N",       100,   "Number of spatial points for output file.");
    opt.add_size_option   ("nst,s",      1,     "Number of states to find.");
    opt.add_char_option   ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, heavy holes or light holes respectively.");
    opt.add_numeric_option("vcb",        0.00,  "Band-edge potential [eV]");
    opt.add_numeric_option("alpha",      0.00,  "Non-parabolicity parameter [eV^{-1}]");
    opt.add_numeric_option("E-cutoff",          "Cut-off energy for solutions [meV]");

    std::string doc("Find the eigenstates of an infinite quantum well.");

    std::string details("The following output text files are created:\n"
                        "  'E*.r'   \tEnergy of each state:\n"
                        "           \tCOLUMN 1: state index.\n"
                        "           \tCOLUMN 2: energy [meV].\n"
                        "  'wf_*i.r'\tWave function amplitude at each position\n"
                        "           \tCOLUMN 1: position [m].\n"
                        "           \tCOLUMN 2: wave function amplitude [m^{-1/2}].\n"
                        "\n"
                        "\tIn each case, the '*' is replaced by the particle ID and the 'i' is replaced by the number of the state.\n"
                        "\n"
                        "Examples:\n"
                        "   Compute the ground state in a 150-angstrom well with effective mass = 0.1 m0:\n\n"
                        "   efiw --width 150 --mass 0.1\n"
                        "\n"
                        "   Compute the first three heavy-hole states in a 200-angstrom well, using effective mass = 0.62 m0:\n\n"
                        "   efiw --width 200 --mass 0.62 --particle h");

    opt.add_prog_specific_options_and_parse(argc, argv, doc, details);

    return opt;
};

int main(int argc, char *argv[])
{
    Options opt = configure_options(argc, argv);

    const double L     = opt.get_numeric_option("width") * 1e-10; // well width [m]
    const char   p     = opt.get_char_option("particle");         // particle ID (e, h or l)
    const double m     = opt.get_numeric_option("mass") * me;     // effective mass [kg]
    const size_t N     = opt.get_size_option("nz");               // number of spatial steps
    const size_t s     = opt.get_size_option("nst");              // number of states
    const double alpha = opt.get_numeric_option("alpha") / e;     // Non-parabolicity [1/J]
    const double V     = opt.get_numeric_option("vcb") * e;       // band_edge potential [J]

    SchroedingerSolverInfWell se(m, L, N, alpha, V, s);

    // Set cut-off energy if desired
    if(opt.vm.count("E-cutoff") > 0)
        se.set_E_cutoff(opt.get_numeric_option("E-cutoff") * e/1000);

    std::vector<State> solutions = se.get_solutions(true);

    // Dump to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);
    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         solutions,
                         se.get_z(),
                         true);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
